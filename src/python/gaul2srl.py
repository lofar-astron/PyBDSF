
"""Module gaul2srl

This will group gaussians in an island into sources. Will code callgaul2srl.f here, though
it could probably be made more efficient.

img.sources is a list of source objects, which are instances of the class Source
(with attributes the same as in .srl of fbdsm).
img.sources[n] is a source.
source.gaussians is the list of component gaussian objects.
source.island_id is the island id of that source.
source.source_id is the source id of that source, the index of source in img.sources.
Each gaussian object gaus has gaus.source_id, the source id.

Also, each island object of img.islands list has the source object island.source
"""

from image import *
from islands import *
from gausfit import Gaussian
from interface import wrap
import mylogger
import numpy as N
N.seterr(divide='raise')

nsrc = Int(doc="Number of sources in the image")
Gaussian.source_id = Int(doc="Source number of a gaussian", colname='Source_id')
Gaussian.code = String(doc='Source code S, C, or M', colname='S_Code')

class Op_gaul2srl(Op):
    """
    Slightly modified from fortran.
    """

    def __call__(self, img):
        #  for each island, get the gaussians into a list and then send them to process
        #  src_index is source number, starting from 0
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Gaul2Srl")
        mylogger.userinfo(mylog, 'Grouping Gaussians into sources')
        img.aperture = img.opts.aperture
        if img.aperture != None and img.aperture <= 0.0:
            mylog.warn('Specified aperture is <= 0. Skipping aperture fluxes.')
            img.aperture = None

        src_index = -1
        dsrc_index = 0
        sources = []
        dsources = []
        no_gaus_islands = []
        for iisl, isl in enumerate(img.islands):
            isl_sources = []
            isl_dsources = []
            g_list = []
            for g in isl.gaul:
                if g.flag == 0:
                    g_list.append(g)

            if len(g_list) > 0:
              if len(g_list) == 1:
                src_index, source = self.process_single_gaussian(img, g_list, src_index, code = 'S')
                sources.append(source)
                isl_sources.append(source)
              else:
                src_index, source = self.process_CM(img, g_list, isl, src_index)
                sources.extend(source)
                isl_sources.extend(source)
            else:
                if not img.waveletimage:
                    dg = isl.dgaul[0]
                    no_gaus_islands.append((isl.island_id, dg.centre_pix[0], dg.centre_pix[1]))
                    # Put in the dummy Source as the source and use negative IDs
                    g_list = isl.dgaul
                    dsrc_index, dsource = self.process_single_gaussian(img, g_list, dsrc_index, code = 'S')
                    dsources.append(dsource)
                    isl_dsources.append(dsource)

            isl.sources = isl_sources
            isl.dsources = isl_dsources
        img.sources = sources
        img.dsources = dsources
        img.nsrc = src_index + 1
        mylogger.userinfo(mylog, "Number of sources formed from Gaussians",
                          str(img.nsrc))
        if not img.waveletimage and not img._pi and len(no_gaus_islands) > 0 and not img.opts.quiet:
            message = 'All Gaussians were flagged for the following island'
            if len(no_gaus_islands) == 1:
                message += ':\n'
            else:
                message += 's:\n'
            for isl_id in no_gaus_islands:
                message += '    Island #%i (x=%i, y=%i)\n' % isl_id
            if len(no_gaus_islands) == 1:
                message += 'Please check this island. If it is a valid island and\n'
            else:
                message += 'Please check these islands. If they are valid islands and\n'
            if img.opts.atrous_do:
                message += 'should be fit, try adjusting the flagging options (use\n'\
                           'show_fit with "ch0_flagged=True" to see the flagged Gaussians).'
            else:
                message += 'should be fit, try adjusting the flagging options (use\n'\
                           'show_fit with "ch0_flagged=True" to see the flagged Gaussians)\n'\
                           'or enabling the wavelet module (with "atrous_do=True").'
            message += '\nTo include these islands in output source catalogs, set\n'\
                        'incl_empty=True in the write_catalog task.'
            mylog.warning(message)

        img.completed_Ops.append('gaul2srl')

#################################################################################################

    def process_single_gaussian(self, img, g_list, src_index, code):
        """ Process single gaussian into a source, for both S and C type sources. g is just one
            Gaussian object (not a list)."""

        g = g_list[0]

        total_flux = [g.total_flux, g.total_fluxE]
        peak_flux_centroid = peak_flux_max = [g.peak_flux, g.peak_fluxE]
        posn_sky_centroid = posn_sky_max = [g.centre_sky, g.centre_skyE]
        size_sky = [g.size_sky, g.size_skyE]
        deconv_size_sky = [g.deconv_size_sky, g.deconv_size_skyE]
        bbox = img.islands[g.island_id].bbox
        ngaus = 1
        island_id = g.island_id
        if g.gaus_num < 0:
            gaussians = []
        else:
            gaussians = list([g])
        aper_flux = func.ch0_aperture_flux(img, g.centre_pix, img.aperture)

        source_prop = list([code, total_flux, peak_flux_centroid, peak_flux_max, aper_flux, posn_sky_centroid, \
             posn_sky_max, size_sky, deconv_size_sky, bbox, ngaus, island_id, gaussians])
        source = Source(img, source_prop)

        if g.gaussian_idx == -1:
            src_index -= 1
        else:
            src_index += 1
        g.source_id = src_index
        g.code = code
        source.source_id = src_index

        return src_index, source

##################################################################################################

    def process_CM(self, img, g_list, isl, src_index):
        """
        Bundle errors with the quantities.
        ngau = number of gaussians in island
        src_id = the source index array for every gaussian in island
        nsrc = final number of distinct sources in the island
        """

        ngau = len(g_list)  # same as cisl in callgaul2srl.f
        nsrc = ngau         # same as islct; initially make each gaussian as a source
        src_id = N.arange(nsrc)  # same as islnum in callgaul2srl.f

        boxx, boxy = isl.bbox
        subn = boxx.stop-boxx.start; subm = boxy.stop-boxy.start
        delc = [boxx.start, boxy.start]
        subim = self.make_subim(subn, subm, g_list, delc)

        index = [(i,j) for i in range(ngau) for j in range(ngau) if j > i]
        for pair in index:
            same_island = self.in_same_island(pair, img, g_list, isl, subim, subn, subm, delc)
            if same_island:
              nsrc -= 1
              mmax, mmin = max(src_id[pair[0]],src_id[pair[1]]), min(src_id[pair[0]],src_id[pair[1]])
              arr = N.where(src_id == mmax)[0]; src_id[arr] = mmin
                            # now reorder src_id so that it is contiguous
        for i in range(ngau):
            ind1 = N.where(src_id==i)[0]
            if len(ind1) == 0:
                arr = N.where(src_id > i)[0]
                if len(arr) > 0:
                  decr =  N.min(src_id[arr])-i
                  for j in arr: src_id[j] -= decr
        nsrc = N.max(src_id)+1
        # now do whats in sub_calc_para_source

        source_list = []
        for isrc in range(nsrc):
          posn = N.where(src_id == isrc)[0]
          g_sublist=[]
          for i in posn:
              g_sublist.append(g_list[i])
          ngau_insrc = len(posn)
                                # Do source type C
          if ngau_insrc == 1:
              src_index, source = self.process_single_gaussian(img, g_sublist, src_index, code = 'C')
          else:
              # make mask and subim. Invalid mask value is -1 since 0 is valid srcid
              mask = self.make_mask(isl, subn, subm, 1, isrc, g_sublist, delc)
              src_index, source = self.process_Multiple(img, g_sublist, mask, src_index, isrc, subim, \
                                  isl, delc, subn, subm)
          source_list.append(source)

        return src_index, source_list

##################################################################################################

    def in_same_island(self, pair, img, g_list, isl, subim, subn, subm, delc):
        """ Whether two gaussians belong to the same source or not. """
        import functions as func

        def same_island_min(pair, g_list, subim, delc, tol=0.5):
            """ If the minimum of the reconstructed fluxes along the line joining the peak positions
                is greater than thresh_isl times the rms_clip, they belong to different islands. """

            g1 = g_list[pair[0]]
            g2 = g_list[pair[1]]
            pix1 = N.array(g1.centre_pix)
            pix2 = N.array(g2.centre_pix)

            x1, y1 = N.floor(pix1)-delc; x2, y2 = N.floor(pix2)-delc
            pix1 = N.array(N.unravel_index(N.argmax(subim[x1:x1+2,y1:y1+2]), (2,2)))+[x1,y1]
            pix2 = N.array(N.unravel_index(N.argmax(subim[x2:x2+2,y2:y2+2]), (2,2)))+[x2,y2]
            if pix1[1] >= subn: pix1[1] = pix1[1]-1
            if pix2[1] >= subm: pix2[1] = pix2[1]-1

            maxline = int(round(N.max(N.abs(pix1-pix2)+1)))
            flux1 = g1.peak_flux
            flux2 = g2.peak_flux
                                                # get pix values of the line
            pixdif = pix2 - pix1
            same_island_min = False
            same_island_cont = False
            if maxline == 1:
              same_island_min = True
              same_island_cont = True
            else:
              if abs(pixdif[0]) > abs(pixdif[1]):
                xline = N.round(min(pix1[0],pix2[0])+N.arange(maxline))
                yline = N.round((pix1[1]-pix2[1])/(pix1[0]-pix2[0])* \
                       (min(pix1[0],pix2[0])+N.arange(maxline)-pix1[0])+pix1[1])
              else:
                yline = N.round(min(pix1[1],pix2[1])+N.arange(maxline))
                xline = N.round((pix1[0]-pix2[0])/(pix1[1]-pix2[1])* \
                       (min(pix1[1],pix2[1])+N.arange(maxline)-pix1[1])+pix1[0])
              rpixval = N.zeros(maxline)
              xbig = N.where(xline >= N.size(subim,0))
              xline[xbig] = N.size(subim,0) - 1
              ybig = N.where(yline >= N.size(subim,1))
              yline[ybig] = N.size(subim,1) - 1
              for i in range(maxline):
                pixval = subim[xline[i],yline[i]]
                rpixval[i] = pixval
              min_pixval = N.min(rpixval)
              minind_p = N.argmin(rpixval)
              maxind_p = N.argmax(rpixval)

              if minind_p in (0, maxline-1) and maxind_p in (0, maxline-1):
                same_island_cont = True
              if min_pixval >= min(flux1, flux2):
                same_island_min = True
              elif abs(min_pixval-min(flux1,flux2)) <= tol*isl.rms*img.opts.thresh_isl:
                same_island_min = True

            return same_island_min, same_island_cont

        def same_island_dist(pair, g_list, tol=0.5):
            """ If the centres are seperated by a distance less than half the sum of their
                fwhms along the PA of the line joining them, they belong to the same island. """
            from math import sqrt

            g1 = g_list[pair[0]]
            g2 = g_list[pair[1]]
            pix1 = N.array(g1.centre_pix)
            pix2 = N.array(g2.centre_pix)
            gsize1 = g1.size_pix
            gsize2 = g2.size_pix

            fwhm1 = func.gdist_pa(pix1, pix2, gsize1)
            fwhm2 = func.gdist_pa(pix1, pix2, gsize2)
            dx = pix2[0]-pix1[0]; dy = pix2[1]-pix1[1]
            dist = sqrt(dy*dy + dx*dx)

            if dist <= tol*(fwhm1+fwhm2):
              same_island = True
            else:
              same_island = False

            return same_island

        if img.opts.group_by_isl:
            same_isl1_min = True
            same_isl1_cont = True
            same_isl2 = True
        else:
            if img.opts.group_method == 'curvature':
                subim = -1.0 * func.make_curvature_map(subim)
            tol = img.opts.group_tol
            same_isl1_min, same_isl1_cont = same_island_min(pair, g_list, subim, delc, tol)
            same_isl2 = same_island_dist(pair, g_list, tol/2.0)

        g1 = g_list[pair[0]]

        same_island = (same_isl1_min and same_isl2) or same_isl1_cont

        return same_island

##################################################################################################

    def process_Multiple(self, img, g_sublist, mask, src_index, isrc, subim, isl, delc, subn, subm):
        """ Same as gaul_to_source.f. isrc is same as k in the fortran version. """
        from math import pi, sqrt
        from const import fwsig
        from scipy import ndimage
        import functions as func

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Gaul2Srl  ")
        dum = img.beam[0]*img.beam[1]
        cdeltsq = img.wcs_obj.acdelt[0]*img.wcs_obj.acdelt[1]
        bmar_p = 2.0*pi*dum/(cdeltsq*fwsig*fwsig)

                                        # try
        subim_src = self.make_subim(subn, subm, g_sublist, delc)
        mompara = func.momanalmask_gaus(subim_src, mask, isrc, bmar_p, True)
                                        # initial peak posn and value
        maxv = N.max(subim_src)
        maxx, maxy = N.unravel_index(N.argmax(subim_src), subim_src.shape)
                                        # fit gaussian around this posn
        blc = N.zeros(2); trc = N.zeros(2)
        n, m = subim_src.shape[0:2]
        bm_pix = N.array([img.pixel_beam()[0]*fwsig, img.pixel_beam()[1]*fwsig, img.pixel_beam()[2]])
        ssubimsize = max(N.int(N.round(N.max(bm_pix[0:2])*2))+1, 5)
        blc[0] = max(0, maxx-(ssubimsize-1)/2); blc[1] = max(0, maxy-(ssubimsize-1)/2)
        trc[0] = min(n, maxx+(ssubimsize-1)/2); trc[1] = min(m, maxy+(ssubimsize-1)/2)
        s_imsize = trc - blc + 1

        p_ini = [maxv, (s_imsize[0]-1)/2.0*1.1, (s_imsize[1]-1)/2.0*1.1, bm_pix[0]/fwsig*1.3, \
                 bm_pix[1]/fwsig*1.1, bm_pix[2]*2]
        data = subim_src[blc[0]:blc[0]+s_imsize[0], blc[1]:blc[1]+s_imsize[1]]
        smask = mask[blc[0]:blc[0]+s_imsize[0], blc[1]:blc[1]+s_imsize[1]]
        rmask = N.where(smask==isrc, False, True)
        x_ax, y_ax = N.indices(data.shape)

        if N.sum(~rmask) >=6:
          para, ierr = func.fit_gaus2d(data, p_ini, x_ax, y_ax, rmask)
          if (0.0<para[1]<s_imsize[0]) and (0.0<para[2]<s_imsize[1]) and \
            para[3]<s_imsize[0] and para[4]<s_imsize[1]:
            maxpeak = para[0]
          else:
            maxpeak = maxv
          posn = para[1:3]-(0.5*N.sum(s_imsize)-1)/2.0+N.array([maxx, maxy])-1+delc
        else:
          maxpeak = maxv
          posn = N.unravel_index(N.argmax(data*~rmask), data.shape)+N.array(delc) +blc

        # calculate peak by bilinear interpolation around centroid
        # First check that moment analysis gave a valid position. If not, use
        # posn from gaussian fit instead.
        if N.isnan(mompara[1]):
            mompara[1] = posn[0] - delc[0]
        x1 = N.int(N.floor(mompara[1]))
        if N.isnan(mompara[2]):
            mompara[2] = posn[1] - delc[1]
        y1 = N.int(N.floor(mompara[2]))
        xind = slice(x1, x1+2, 1); yind = slice(y1, y1+2, 1)
        if img.opts.flag_smallsrc and (N.sum(mask[xind, yind]==N.ones((2,2))*isrc) != 4):
            mylog.debug('Island = '+str(isl.island_id))
            mylog.debug('Mask = '+repr(mask[xind, yind])+'xind, yind, x1, y1 = '+repr(xind)+' '+repr(yind)+' '+repr(x1)+' '+repr(y1))
        t=(mompara[1]-x1)/(x1+1-x1)  # in case u change it later
        u=(mompara[2]-y1)/(y1+1-y1)
        s_peak=(1.0-t)*(1.0-u)*subim_src[x1,y1]+t*(1.0-u)*subim_src[x1+1,y1]+ \
               t*u*subim_src[x1+1,y1+1]+(1.0-t)*u*subim_src[x1,y1+1]
        if (not img.opts.flag_smallsrc) and (N.sum(mask[xind, yind]==N.ones((2,2))*isrc) != 4):
            mylog.debug('Speak '+repr(s_peak)+'Mompara = '+repr(mompara))
            mylog.debug('x1, y1 : '+repr(x1)+', '+repr(y1))
            # import pylab as pl
            # pl.imshow(N.transpose(subim_src), origin='lower', interpolation='nearest')
            # pl.suptitle('Image of bad M source '+str(isl.island_id))
                                        # convert pixels to coords
        try:
            sra, sdec = img.pix2sky([mompara[1]+delc[0], mompara[2]+delc[1]])
            mra, mdec = img.pix2sky(posn)
        except RuntimeError, err:
            # Invalid pixel wcs coordinate
            sra, sdec = 0.0, 0.0
            mra, mdec = 0.0, 0.0
                                        # "deconvolve" the sizes
        gaus_c = [mompara[3], mompara[4], mompara[5]]
        gaus_bm = [bm_pix[0], bm_pix[1], bm_pix[2]]
        gaus_dc, err = func.deconv2(gaus_bm, gaus_c)
        deconv_size_sky = img.pix2gaus(gaus_dc, [mompara[1]+delc[0], mompara[2]+delc[1]])

                                        # update all objects etc
        tot = 0.0
        totE_sq = 0.0
        for g in g_sublist:
            tot += g.total_flux
            totE_sq += g.total_fluxE**2
        totE = sqrt(totE_sq)
        size_pix = [mompara[3], mompara[4], mompara[5]]
        size_sky = img.pix2gaus(size_pix, [mompara[1]+delc[0], mompara[2]+delc[1]])

        # Estimate uncertainties in source size and position due to
        # errors in the constituent Gaussians using a Monte Carlo technique.
        # Sum with Condon (1997) errors in quadrature.
        plist = mompara.tolist()+[tot]
        plist[0] = s_peak
        plist[3] /= fwsig
        plist[4] /= fwsig
        errors = func.get_errors(img, plist, isl.rms)

        if img.opts.do_mc_errors:
            nMC = 20
            mompara0_MC = N.zeros(nMC, dtype=float)
            mompara1_MC = N.zeros(nMC, dtype=float)
            mompara2_MC = N.zeros(nMC, dtype=float)
            mompara3_MC = N.zeros(nMC, dtype=float)
            mompara4_MC = N.zeros(nMC, dtype=float)
            mompara5_MC = N.zeros(nMC, dtype=float)
            for i in range(nMC):
                # Reconstruct source from component Gaussians. Draw the Gaussian
                # parameters from random distributions given by their errors.
                subim_src_MC = self.make_subim(subn, subm, g_sublist, delc, mc=True)

                try:
                    mompara_MC = func.momanalmask_gaus(subim_src_MC, mask, isrc, bmar_p, True)
                    mompara0_MC[i] = mompara_MC[0]
                    mompara1_MC[i] = mompara_MC[1]
                    mompara2_MC[i] = mompara_MC[2]
                    mompara3_MC[i] = mompara_MC[3]
                    mompara4_MC[i] = mompara_MC[4]
                    mompara5_MC[i] = mompara_MC[5]
                except:
                    mompara0_MC[i] = mompara[0]
                    mompara1_MC[i] = mompara[1]
                    mompara2_MC[i] = mompara[2]
                    mompara3_MC[i] = mompara[3]
                    mompara4_MC[i] = mompara[4]
                    mompara5_MC[i] = mompara[5]
            mompara0E = N.std(mompara0_MC)
            mompara1E = N.std(mompara1_MC)
            if mompara1E > 2.0*mompara[1]:
                mompara1E = 2.0*mompara[1] # Don't let errors get too large
            mompara2E = N.std(mompara2_MC)
            if mompara2E > 2.0*mompara[2]:
                mompara2E = 2.0*mompara[2] # Don't let errors get too large
            mompara3E = N.std(mompara3_MC)
            if mompara3E > 2.0*mompara[3]:
                mompara3E = 2.0*mompara[3] # Don't let errors get too large
            mompara4E = N.std(mompara4_MC)
            if mompara4E > 2.0*mompara[4]:
                mompara4E = 2.0*mompara[4] # Don't let errors get too large
            mompara5E = N.std(mompara5_MC)
            if mompara5E > 2.0*mompara[5]:
                mompara5E = 2.0*mompara[5] # Don't let errors get too large
        else:
             mompara1E = 0.0
             mompara2E = 0.0
             mompara3E = 0.0
             mompara4E = 0.0
             mompara5E = 0.0

        # Now add MC errors in quadrature with Condon (1997) errors
        size_skyE = [sqrt(mompara3E**2 + errors[3]**2) * sqrt(cdeltsq),
                     sqrt(mompara4E**2 + errors[4]**2) * sqrt(cdeltsq),
                     sqrt(mompara5E**2 + errors[5]**2)]
        sraE, sdecE = (sqrt(mompara1E**2 + errors[1]**2) * sqrt(cdeltsq),
                       sqrt(mompara2E**2 + errors[2]**2) * sqrt(cdeltsq))
        deconv_size_skyE = size_skyE # set deconvolved errors to non-deconvolved ones

        # Find aperture flux
        if img.opts.aperture_posn == 'centroid':
            aper_pos = [mompara[1]+delc[0], mompara[2]+delc[1]]
        else:
            aper_pos = posn
        aper_flux, aper_fluxE = func.ch0_aperture_flux(img, aper_pos, img.aperture)

        isl_id = isl.island_id
        source_prop = list(['M', [tot, totE], [s_peak, isl.rms], [maxpeak, isl.rms],
                      [aper_flux, aper_fluxE], [[sra, sdec],
                      [sraE, sdecE]], [[mra, mdec], [sraE, sdecE]], [size_sky, size_skyE],
                      [deconv_size_sky, deconv_size_skyE], isl.bbox, len(g_sublist),
                      isl_id, g_sublist])
        source = Source(img, source_prop)

        src_index += 1
        for g in g_sublist:
            g.source_id = src_index
            g.code = 'M'
        source.source_id = src_index

        return src_index, source

##################################################################################################

    def make_subim(self, subn, subm, g_list, delc, mc=False):
        import functions as func

        subim = N.zeros((subn, subm))
        x, y = N.indices((subn, subm))
        for g in g_list:
            params = func.g2param(g)
            params[1] -= delc[0]; params[2] -= delc[1]
            if mc:
                # draw random variables from distributions given by errors
                params_err = func.g2param_err(g)
                for i in range(len(params)):
                    mc_param = N.random.normal(loc=params[i], scale=params_err[i])
                    params[i] = mc_param
            gau = func.gaus_2d(params, x, y)
            subim = subim + gau

        return subim

##################################################################################################

    def make_mask(self, isl, subn, subm, nsrc, src_id, g_list, delc):
        import functions as func
                                        # define stuff for calculating gaussian
        boxx, boxy = isl.bbox
        subn = boxx.stop-boxx.start; subm = boxy.stop-boxy.start
        x, y = N.indices((subn, subm))
                                        # construct image of each source in the island
        src_image = N.zeros((subn, subm, nsrc))
        nn = 1
        for isrc in range(nsrc):
            if nsrc == 1:
                g_sublist = g_list
            else:
                posn = N.where(src_id == isrc)[0]
                g_sublist=[]
                for i in posn:
                    g_sublist.append(g_list[i])
            for g in g_sublist:
                params = func.g2param(g)
                params[1] -= delc[0]; params[2] -= delc[1]
                gau = func.gaus_2d(params, x, y)
                src_image[:,:,isrc] = src_image[:,:,isrc] + gau
                                        # mark each pixel as belonging to one source
                                        # just compare value, should compare with sigma later
        mask = N.argmax(src_image, axis=2) + src_id
        orig_mask = isl.mask_active
        mask[N.where(orig_mask)] = -1

        return mask


##################################################################################################
#  Define class Source
##################################################################################################

from image import *
from gausfit import Gaussian
from islands import Island

class Source(object):
    """ Instances of this class store sources made from grouped gaussians. """
    source_id           = Int(doc="Source index", colname='Source_id')
    code                = String(doc='Source code S, C, or M', colname='S_Code')
    total_flux          = Float(doc="Total flux density (Jy)", colname='Total_flux', units='Jy')
    total_fluxE         = Float(doc="Error in total flux density (Jy)", colname='E_Total_flux',
                                units='Jy')
    peak_flux_centroid  = Float(doc="Peak flux density per beam at centroid of emission (Jy/beam)",
                                colname='Peak_flux_cen', units='Jy/beam')
    peak_flux_centroidE = Float(doc="Error in peak flux density per beam at centroid of emission (Jy/beam)",
                                colname='E_Peak_flux_cen', units='Jy/beam')
    peak_flux_max       = Float(doc="Peak flux density per beam at posn of maximum emission (Jy/beam)",
                                colname='Peak_flux', units='Jy/beam')
    peak_flux_maxE      = Float(doc="Error in peak flux density per beam at posn of max emission (Jy/beam)",
                                colname='E_Peak_flux', units='Jy/beam')
    aperture_flux       = Float(doc="Total aperture flux density (Jy)", colname='Aperture_flux',
                                units='Jy')
    aperture_fluxE      = Float(doc="Error in total aperture flux density (Jy)", colname='E_Aperture_flux',
                                units='Jy')
    posn_sky_centroid   = List(Float(), doc="Posn (RA, Dec in deg) of centroid of source",
                               colname=['RA', 'DEC'], units=['deg', 'deg'])
    posn_sky_centroidE  = List(Float(), doc="Error in posn (RA, Dec in deg) of centroid of source",
                               colname=['E_RA', 'E_DEC'], units=['deg', 'deg'])
    posn_sky_max        = List(Float(), doc="Posn (RA, Dec in deg) of maximum emission of source",
                               colname=['RA_max', 'DEC_max'], units=['deg', 'deg'])
    posn_sky_maxE       = List(Float(), doc="Error in posn (deg) of maximum emission of source",
                               colname=['E_RA_max', 'E_DEC_max'], units=['deg', 'deg'])
    posn_pix_centroid   = List(Float(), doc="Position (x, y in pixels) of centroid of source",
                               colname=['Xposn', 'Yposn'], units=['pix', 'pix'])
    posn_pix_centroidE  = List(Float(), doc="Error in position (x, y in pixels) of centroid of source",
                               colname=['E_Xposn', 'E_Yposn'], units=['pix', 'pix'])
    posn_pix_max        = List(Float(), doc="Position (x, y in pixels) of maximum emission of source",
                               colname=['Xposn_max', 'Yposn_max'], units=['pix', 'pix'])
    posn_pix_maxE       = List(Float(), doc="Error in position (pixels) of maximum emission of source",
                               colname=['E_Xposn_max', 'E_Yposn_max'], units=['pix', 'pix'])
    size_sky            = List(Float(), doc="Shape of the source FWHM, BPA, deg",
                               colname=['Maj', 'Min', 'PA'], units=['deg', 'deg',
                              'deg'])
    size_skyE           = List(Float(), doc="Error on shape of the source FWHM, BPA, deg",
                               colname=['E_Maj', 'E_Min', 'E_PA'], units=['deg', 'deg',
                               'deg'])
    deconv_size_sky     = List(Float(), doc="Deconvolved shape of the source FWHM, BPA, deg",
                               colname=['DC_Maj', 'DC_Min', 'DC_PA'], units=['deg', 'deg',
                              'deg'])
    deconv_size_skyE    = List(Float(), doc="Error on deconvolved shape of the source FWHM, BPA, deg",
                               colname=['E_DC_Maj', 'E_DC_Min', 'E_DC_PA'], units=['deg', 'deg',
                              'deg'])
    rms_isl             = Float(doc="Island rms Jy/beam", colname='Isl_rms', units='Jy/beam')
    mean_isl            = Float(doc="Island mean Jy/beam", colname='Isl_mean', units='Jy/beam')
    total_flux_isl      = Float(doc="Island total flux from sum of pixels", colname='Isl_Total_flux', units='Jy')
    total_flux_islE     = Float(doc="Error on island total flux from sum of pixels", colname='E_Isl_Total_flux', units='Jy')
    gresid_rms          = Float(doc="Island rms in Gaussian residual image Jy/beam",
                                colname='Resid_Isl_rms', units='Jy/beam')
    gresid_mean         = Float(doc="Island mean in Gaussian residual image Jy/beam",
                                colname='Resid_Isl_mean', units='Jy/beam')
    sresid_rms          = Float(doc="Island rms in Shapelet residual image Jy/beam",
                                colname='Resid_Isl_rms', units='Jy/beam')
    sresid_mean         = Float(doc="Island mean in Shapelet residual image Jy/beam",
                                colname='Resid_Isl_mean', units='Jy/beam')
    ngaus               = Int(doc='Number of gaussians in the source', colname='N_gaus')
    island_id           = Int(doc="Serial number of the island", colname='Isl_id')
    gaussians           = List(tInstance(Gaussian), doc="")
    bbox                = List(Instance(slice(0), or_none=False), doc = "")

    def __init__(self, img, sourceprop):

        code, total_flux, peak_flux_centroid, peak_flux_max, aper_flux, posn_sky_centroid, \
                     posn_sky_max, size_sky, deconv_size_sky, bbox, ngaus, island_id, gaussians = sourceprop
        self.code = code
        self.total_flux, self.total_fluxE = total_flux
        self.peak_flux_centroid, self.peak_flux_centroidE = peak_flux_centroid
        self.peak_flux_max, self.peak_flux_maxE = peak_flux_max
        self.posn_sky_centroid, self.posn_sky_centroidE = posn_sky_centroid
        self.posn_sky_max, self.posn_sky_maxE = posn_sky_max
        self.size_sky, self.size_skyE = size_sky
        self.deconv_size_sky, self.deconv_size_skyE = deconv_size_sky
        self.bbox = bbox
        self.ngaus = ngaus
        self.island_id = island_id
        self.gaussians = gaussians
        self.rms_isl = img.islands[island_id].rms
        self.mean_isl = img.islands[island_id].mean
        self.total_flux_isl = img.islands[island_id].total_flux
        self.total_flux_islE = img.islands[island_id].total_fluxE
        self.mean_isl = img.islands[island_id].mean
        self.jlevel = img.j
        self.aperture_flux, self.aperture_fluxE =  aper_flux


Image.sources = List(tInstance(Source), doc="List of Sources")
Island.sources = List(tInstance(Source), doc="List of Sources")



