"""Module gausfit.

This module does multi-gaussian fits for all detected islands.
At the moment fitting algorithm is quite simple -- we just add
gaussians one-by-one as long as there are pixels with emission
in the image, and do post-fitting flagging of the extracted
gaussians.

The fitting itself is implemented by the means of MGFunction
class and a number of fitter routines in _cbdsm module.
MGFunction class implements multi-gaussian function and
provides all functionality required by the specific fitters.
"""

from image import *
from copy import deepcopy as cp
import mylogger
import sys
import time
import statusbar
import _cbdsm
try:
    import matplotlib.pyplot as pl
    has_pl = True
except ImportError:
    has_pl = False
import scipy.ndimage as nd
import multi_proc as mp


ngaus = Int(doc="Total number of gaussians extracted")
total_flux_gaus = Float(doc="Total flux in the Gaussians extracted")

class Op_gausfit(Op):
    """Fit a number of 2D gaussians to each island.

    The results of the fitting are stored in the Island
    structure itself as a list of Gaussian objects (gaul) and a
    list of flagged gaussians (fgaul).

    Prerequisites: module islands should be run first.
    """
    def __call__(self, img):
        from time import time
        import functions as func
        import itertools

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Gausfit")
        if len(img.islands) == 0:
            img.gaussians = []
            img.ngaus = 0
            img.total_flux_gaus = 0.0
            img.completed_Ops.append('gausfit')
            return img

        bar = statusbar.StatusBar('Fitting islands with Gaussians .......... : ',
                                  0, img.nisl)
        opts = img.opts
        if opts.quiet == False and opts.verbose_fitting == False:
            bar.start()
        iter_ngmax  = 10
        min_maxsize = 50.0
        maxsize = opts.splitisl_maxsize
        min_peak_size = 30.0
        peak_size = opts.peak_maxsize
        if maxsize < min_maxsize:
            maxsize = min_maxsize
            opts.splitisl_maxsize = min_maxsize
        if peak_size < min_peak_size:
            peak_size = min_peak_size
            opts.peak_maxsize = min_peak_size

        # Set up multiproccessing. First create a simple copy of the Image
        # object that contains the minimal data needed.
        opts_dict = opts.to_dict()
        img_simple = Image(opts_dict)
        img_simple.pixel_beamarea = img.pixel_beamarea
        img_simple.pixel_beam = img.pixel_beam
        img_simple.thresh_pix = img.thresh_pix
        img_simple.minpix_isl = img.minpix_isl
        img_simple.clipped_mean = img.clipped_mean

        # Next, define the weights to use when distributing islands among cores.
        # The weight should scale with the processing time. At the moment
        # we use the island area, but other parameters may be better.
        weights = []
        for isl in img.islands:
            weights.append(isl.size_active)

        # Now call the parallel mapping function. Returns a list of [gaul, fgaul]
        # for each island.
        gaus_list = mp.parallel_map(func.eval_func_tuple,
                    itertools.izip(itertools.repeat(self.process_island),
                    img.islands, itertools.repeat(img_simple),
                    itertools.repeat(opts)), numcores=opts.ncores,
                    bar=bar, weights=weights)

        for isl in img.islands:
            ### now convert gaussians into Gaussian objects and store
            idx = isl.island_id
            gaul = gaus_list[idx][0]
            fgaul = gaus_list[idx][1]
            gaul = [Gaussian(img, par, idx, gidx)
                        for (gidx, par) in enumerate(gaul)]

            if len(gaul) == 0:
                gidx = 0
            fgaul= [Gaussian(img, par, idx, gidx + gidx2 + 1, flag)
                        for (gidx2, (flag, par)) in enumerate(fgaul)]

            isl.gaul = gaul
            isl.fgaul= fgaul

        gaussian_list = [g for isl in img.islands for g in isl.gaul]
        img.gaussians = gaussian_list

        ### put in the serial number of the gaussians for the whole image
        n = 0
        tot_flux = 0.0
        if img.waveletimage:
            # store the wavelet scale for each Gaussian
            # (wavelet img's have a img.j attribute)
            j = img.j
        else:
            j = 0
        for isl in img.islands:
            m = 0
            for g in isl.gaul:
                n += 1; m += 1
                g.gaus_num = n - 1
                tot_flux += g.total_flux
            isl.ngaus = m
        img.ngaus = n
        img.total_flux_gaus = tot_flux

        mylogger.userinfo(mylog, "Total number of Gaussians fit to image",
                          str(n))
        if not hasattr(img, '_pi') and not img.waveletimage:
            mylogger.userinfo(mylog, "Total flux density in model", '%.3f Jy' %
                          tot_flux)

        # Check if model flux is very different from sum of flux in image
        if img.ch0_sum_jy > 0 and not hasattr(img, '_pi'):
            if img.total_flux_gaus/img.ch0_sum_jy < 0.5 or \
                    img.total_flux_gaus/img.ch0_sum_jy > 2.0:
                mylog.warn('Total flux density in model is %0.2f times sum of pixels '\
                               'in input image. Large residuals may remain.' %
                           (img.total_flux_gaus/img.ch0_sum_jy,))

        # Check if there are many Gaussians with deconvolved size of 0 in one
        # axis but not in the other. Don't bother to do this for wavelet images.
        fraction_1d = self.check_for_1d_gaussians(img)
        if fraction_1d > 0.5 and img.beam != None and img.waveletimage == False:
            mylog.warn('After deconvolution, more than 50% of Gaussians are '\
                           "1-D. Unless you're fitting an extended source, "\
                           "beam may be incorrect.")

        img.completed_Ops.append('gausfit')
        return img


    def process_island(self, isl, img, opts=None):
        """Processes a single island.

        Returns a list best-fit Gaussians and flagged Gaussians.
        """
        import functions as func

        if opts == None:
            opts = img.opts
        iter_ngmax  = 10
        maxsize = opts.splitisl_maxsize
        min_peak_size = 30.0
        min_maxsize = 50.0
        peak_size = opts.peak_maxsize
        if maxsize < min_maxsize:
            maxsize = min_maxsize
            opts.splitisl_maxsize = min_maxsize
        if peak_size < min_peak_size:
            peak_size = min_peak_size
            opts.peak_maxsize = min_peak_size

        size = isl.size_active/img.pixel_beamarea*2.0   # 2.0 roughly corrects for thresh_isl
        if opts.verbose_fitting:
            print "Fitting isl #", isl.island_id, '; # pix = ',N.sum(~isl.mask_active),'; size = ',size

        if size > maxsize:
            tosplit = func.isl_tosplit(isl, opts)
            if opts.split_isl and tosplit[0] > 0:
                n_subisl, sub_labels = tosplit[1], tosplit[2]
                gaul = []; fgaul = []
                if opts.verbose_fitting:
                    print 'SPLITTING ISLAND INTO ',n_subisl,' PARTS FOR ISLAND ',isl.island_id
                for i_sub in range(n_subisl):
                    islcp = isl.copy(img.pixel_beamarea)
                    islcp.mask_active = N.where(sub_labels == i_sub+1, False, True)
                    islcp.mask_noisy = N.where(sub_labels == i_sub+1, False, True)
                    size_subisl = (~islcp.mask_active).sum()/img.pixel_beamarea*2.0
                    if opts.peak_fit and size_subisl > peak_size:
                        sgaul, sfgaul = self.fit_island_iteratively(img, islcp, iter_ngmax=iter_ngmax, opts=opts)
                    else:
                        sgaul, sfgaul = self.fit_island(islcp, opts, img)
                    gaul = gaul + sgaul; fgaul = fgaul + sfgaul
            else:
                isl.islmean = 0.0
                if opts.peak_fit and size > peak_size:
                    gaul, fgaul = self.fit_island_iteratively(img, isl, iter_ngmax=iter_ngmax, opts=opts)
                else:
                    gaul, fgaul = self.fit_island(isl, opts, img)

        else:
            if opts.peak_fit and size > peak_size:
                gaul, fgaul = self.fit_island_iteratively(img, isl, iter_ngmax=iter_ngmax, opts=opts)
            else:
                gaul, fgaul = self.fit_island(isl, opts, img)

        # Return list of Gaussians
        return [gaul, fgaul]

    def fit_island(self, isl, opts, img, ngmax=None, ffimg=None, ini_gausfit=None):
        """Fit island with a set of 2D gaussians.

        Parameters:
        isl: island
        opts: Opts structure of the image
        beam: beam parameters which are used as an initial guess for
              gaussian shape

        Returns:
        Function returns 2 lists with parameters of good and flagged
        gaussians. Gaussian parameters are updated to be image-relative.

        My own notes (Niruj)
        fcn = MGFunction(im, mask, 1) makes an fcn object
        fcn.find_peak() finds peak and posn in im after subtracting all gaussians in fcn.parameters
        fcn.add_gaussian(gtype, (blah)) adds to fcn.parameters.
        fit(fcn, 0, 1) fits using fcn.parameters as initial guess and overwrites to fcn.parameters.
        fcn.reset() resets just the fcn.parameters. Image is still there.
        Atleast, thats what I think it is.

        """
        from _cbdsm import MGFunction
        if ffimg == None:
            fcn = MGFunction(isl.image-isl.islmean, isl.mask_active, 1)
        else:
            fcn = MGFunction(isl.image-isl.islmean-ffimg, isl.mask_active, 1)
        beam = img.pixel_beam

        if abs(beam[0]/beam[1]) < 1.1:
            beam = (1.1*beam[0], beam[1], beam[2])

        thr1 = isl.mean + opts.thresh_isl*isl.rms
        thr2 = isl.mean + img.thresh_pix*isl.rms
        thr0 = thr1
        verbose = opts.verbose_fitting
        peak = fcn.find_peak()[0]
        dof = isl.size_active
        shape = isl.shape
        isl_image = isl.image - isl.islmean
        size = isl.size_active/img.pixel_beamarea*2.0
        gaul = []
        iter = 0
        ng1 = 0
        if ini_gausfit == None:
            ini_gausfit = opts.ini_gausfit

        if ini_gausfit not in ['default', 'simple', 'nobeam']:
            ini_gausfit = 'default'
        if ini_gausfit == 'simple' and ngmax == None:
          ngmax = 25
        if ini_gausfit == 'default':
          gaul, ng1, ngmax = self.inigaus_fbdsm(isl, thr0, beam, img)
        if ini_gausfit == 'nobeam':
          gaul = self.inigaus_nobeam(isl, thr0, beam, img)
          ng1 = len(gaul); ngmax = ng1+2
        while iter < 5:
            iter += 1
            fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, ini_gausfit, ngmax, verbose)
            gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                              beam, thr0, peak, shape, isl.mask_active,
                                              isl.image, size)
            ng1 = len(gaul)
            if fitok and len(fgaul) == 0:
                break
        if not fitok and ini_gausfit != 'simple':
            # If fits using default or nobeam methods did not work,
            # try using simple instead
            gaul = []
            iter = 0
            ng1 = 0
            ngmax = 25
            while iter < 5:
               iter += 1
               fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, 'simple', ngmax, verbose)
               gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                                 beam, thr0, peak, shape, isl.mask_active,
                                                 isl.image, size)
               ng1 = len(gaul)
               if fitok and len(fgaul) == 0:
                   break
        if not fitok:
            # If normal fitting fails, try to fit 5 or fewer Gaussians to the island
            ngmax = 5
            while not fitok and ngmax > 1:
                fitok = self.fit_iter([], 0, fcn, dof, beam, thr0, 1, 'simple', ngmax, verbose)
                ngmax -= 1
                gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                          beam, thr0, peak, shape, isl.mask_active,
                                          isl.image, size)
        sm_isl = nd.binary_dilation(isl.mask_active)
        if not fitok and N.sum(~sm_isl) >= img.minpix_isl:
            # If all else fails, shrink the island a little and try one last time
            if ffimg == None:
                fcn = MGFunction(isl.image-isl.islmean, nd.binary_dilation(isl.mask_active), 1)
            else:
                fcn = MGFunction(isl.image-isl.islmean-ffimg, nd.binary_dilation(isl.mask_active), 1)
            gaul = []
            iter = 0
            ng1 = 0
            ngmax = 25
            while iter < 5:
               iter += 1
               fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, 'simple', ngmax, verbose)
               gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                                 beam, thr0, peak, shape, isl.mask_active,
                                                 isl.image, size)
               ng1 = len(gaul)
               if fitok and len(fgaul) == 0:
                   break
        lg_isl = nd.binary_erosion(isl.mask_active)
        if not fitok and N.sum(~lg_isl) >= img.minpix_isl:
            # If all else fails, expand the island a little and try one last time
            if ffimg == None:
                fcn = MGFunction(isl.image-isl.islmean, nd.binary_erosion(isl.mask_active), 1)
            else:
                fcn = MGFunction(isl.image-isl.islmean-ffimg, nd.binary_erosion(isl.mask_active), 1)
            gaul = []
            iter = 0
            ng1 = 0
            ngmax = 25
            while iter < 5:
               iter += 1
               fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, 'simple', ngmax, verbose)
               gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                                 beam, thr0, peak, shape, isl.mask_active,
                                                 isl.image, size)
               ng1 = len(gaul)
               if fitok and len(fgaul) == 0:
                   break


        ### return whatever we got
        isl.mg_fcn = fcn
        gaul  = [self.fixup_gaussian(isl, g) for g in gaul]
        fgaul = [(flag, self.fixup_gaussian(isl, g))
                                       for flag, g in fgaul]

        if verbose:
            print 'Number of good Gaussians: %i' % (len(gaul),)
            print 'Number of flagged Gaussians: %i' % (len(fgaul),)
        return gaul, fgaul

    def deblend_and_fit(self, img, isl, opts=None):
        """Deblends an island and then fits it"""
        import functions as func
        sgaul = []; sfgaul = []
        gaul = []; fgaul = []
        if opts == None:
            opts = img.opts
        thresh_isl = opts.thresh_isl
        thresh_pix = opts.thresh_pix
        thresh = opts.fittedimage_clip
        rms = isl.rms
        factor = 1.0
        # Set desired size of sub-island. Don't let it get too small, or fitting
        # won't work well.
        maxsize = max(opts.peak_maxsize, (~isl.mask_active).sum()/img.pixel_beamarea/2.0)

        if opts.verbose_fitting:
            print 'Finding and fitting peaks of island ', isl.island_id
        while True:
            factor *= 1.2
            if N.max(isl.image-isl.islmean-isl.mean)/thresh_isl/factor <= rms:
                if int(factor) == 1:
                    slices = []
                break
            mask_active_orig = isl.mask_active
            act_pixels = (isl.image-isl.islmean-isl.mean)/thresh_isl/factor >= rms
            N.logical_and(act_pixels, ~mask_active_orig, act_pixels)
            rank = len(isl.shape)
            # generates matrix for connectivity, in this case, 8-conn
            connectivity = nd.generate_binary_structure(rank, rank)
            # labels = matrix with value = (initial) island number
            sub_labels, count = nd.label(act_pixels, connectivity)
            # slices has limits of bounding box of each such island
            slices = nd.find_objects(sub_labels)
            if len(slices) == 0:
                break
            size = []
            for idx, s in enumerate(slices):
                idx += 1 # nd.labels indices are counted from 1
                size.append((sub_labels[s] == idx).sum()*2.0)
            # Check whether we have reduced the size of smallest island to less
            # than maxsize; if not, continue with higher threshhold
            if min(size) < maxsize*img.pixel_beamarea:
                break
        gaul = []; fgaul = []
        n_subisl = len(slices)
        if opts.verbose_fitting and n_subisl > 1:
          print 'SEPARATED ISLAND INTO ',n_subisl,' PEAKS FOR ISLAND ',isl.island_id
        for i_sub_isl in range(n_subisl):
          islcp = isl.copy(img)
          islcp.mask_active = N.where(sub_labels == i_sub_isl+1, False, True)
          islcp.mask_noisy = N.where(sub_labels == i_sub_isl+1, False, True)
          sgaul, sfgaul = self.fit_island(islcp, opts, img)
          gaul = gaul + sgaul; fgaul = fgaul + sfgaul
#           if bar.started: bar.spin()

        # Now fit residuals
        ffimg_tot = N.zeros(isl.shape)
        if len(gaul) > 0:
            gaul_obj_list = [Gaussian(img, par, isl.island_id, gidx) for (gidx, par) in enumerate(gaul)]
            for g in gaul_obj_list:
                g.centre_pix[0] -= isl.origin[0]
                g.centre_pix[1] -= isl.origin[1]
                C1, C2 = g.centre_pix
                shape = isl.shape
                b = find_bbox(thresh*isl.rms, g)
                bbox = N.s_[max(0, int(C1-b)):min(shape[0], int(C1+b+1)),
                            max(0, int(C2-b)):min(shape[1], int(C2+b+1))]
                x_ax, y_ax = N.mgrid[bbox]
                ffimg = func.gaussian_fcn(g, x_ax, y_ax)
                ffimg_tot[bbox] += ffimg
        if N.max(isl.image-ffimg_tot-isl.islmean-isl.mean)/thresh_pix >= rms:
            sgaul, sfgaul = self.fit_island(isl, opts, img, ffimg=ffimg_tot)
            gaul = gaul + sgaul; fgaul = fgaul + sfgaul

        return gaul, fgaul

    def fit_island_iteratively(self, img, isl, iter_ngmax=5, opts=None):
        """Fits an island iteratively.

        For large islands, which can require many Gaussians to fit well,
        it is much faster to fit a small number of Gaussians simultaneously
        and iterate."""
        import functions as func
        sgaul = []; sfgaul = []
        gaul = []; fgaul = []
        beam = img.pixel_beam
        if opts == None:
            opts = img.opts
        thresh_isl = opts.thresh_isl
        thresh_pix = opts.thresh_pix
        thresh = opts.fittedimage_clip
        thr = isl.mean + thresh_isl * isl.rms
        rms = isl.rms

        if opts.verbose_fitting:
            print 'Iteratively fitting island ', isl.island_id
        gaul = []; fgaul = []
        ffimg_tot = N.zeros(isl.shape)
        peak_val = N.max(isl.image - isl.islmean)
        while peak_val >= thr:
            sgaul, sfgaul = self.fit_island(isl, opts, img, ffimg=ffimg_tot, ngmax=iter_ngmax, ini_gausfit='simple')
            gaul = gaul + sgaul; fgaul = fgaul + sfgaul

            # Calculate residual image
            if len(sgaul) > 0:
                for g in sgaul:
                    gcopy = g[:]
                    gcopy[1] -= isl.origin[0]
                    gcopy[2] -= isl.origin[1]
                    S1, S2, Th = func.corrected_size(gcopy[3:6])
                    gcopy[3] = S1
                    gcopy[4] = S2
                    gcopy[5] = Th
                    A, C1, C2, S1, S2, Th = gcopy
                    shape = isl.shape
                    b = find_bbox(thresh*isl.rms, gcopy)
                    bbox = N.s_[max(0, int(C1-b)):min(shape[0], int(C1+b+1)),
                                max(0, int(C2-b)):min(shape[1], int(C2+b+1))]
                    x_ax, y_ax = N.mgrid[bbox]
                    ffimg = func.gaussian_fcn(gcopy, x_ax, y_ax)
                    ffimg_tot[bbox] += ffimg
                peak_val = N.max(isl.image - isl.islmean - ffimg_tot)
            else:
                break

        if len(gaul) == 0:
            # Fitting iteratively did not work -- try normal fit
            gaul, fgaul = self.fit_island(isl, opts, img, ini_gausfit='default')

        return gaul, fgaul


    def inigaus_fbdsm(self, isl, thr, beam, img):
        """ initial guess for gaussians like in fbdsm """
        from math import sqrt
        from const import fwsig
        import functions as func

        im = isl.image-isl.islmean; mask = isl.mask_active; av = img.clipped_mean
        inipeak, iniposn, im1 = func.get_maxima(im, mask, thr, isl.shape, beam)
        if len(inipeak) == 0:
          av, stdnew, maxv, maxp, minv, minp = func.arrstatmask(im, mask)
          inipeak = [maxv]; iniposn = [maxp]
        nmulsrc1 = len(iniposn)

        domore = True
        while domore:
          domore = False
          av, stdnew, maxv, maxp, minv, minp = func.arrstatmask(im1, mask)
          if stdnew > isl.rms and maxv >= thr and maxv >= isl.mean+2.0*isl.rms:
            domore = True
            x1, y1 = N.array(iniposn).transpose()
            dumr = N.sqrt((maxp[0]-x1)*(maxp[0]-x1)+(maxp[1]-y1)*(maxp[1]-y1))
            distbm = dumr/sqrt(beam[0]*beam[1]*fwsig*fwsig)
            if N.any((distbm < 0.5) + (dumr < 2.2)):
              domore = False
            if domore:
              iniposn.append(N.array(maxp)); inipeak.append(maxv)
              im1 = func.mclean(im1, maxp, beam)

        inipeak = N.array(inipeak); iniposn = N.array(iniposn)
        ind = list(N.argsort(inipeak)); ind.reverse()
        inipeak = inipeak[ind]
        iniposn = iniposn[ind]
        gaul = []
        for i in range(len(inipeak)):
          g = (float(inipeak[i]), int(iniposn[i][0]), int(iniposn[i][1])) + beam
          gaul.append(g)

        return gaul, nmulsrc1, len(inipeak)

    def inigaus_nobeam(self, isl, thr, beam, img):
        """ To get initial guesses when the source sizes are very different
        from the beam, and can also be elongated. Mainly in the context of
        a-trous transform images. Need to arrive at a good guess of the sizes
        and hence need to partition the image around the maxima first. Tried the
        IFT watershed algo but with markers, it segments the island only around
        the minima and not the whole island. Cant find a good weighting scheme
        for tesselation either. Hence will try this :

        Calculate number of maxima. If one, then take moment as initial
        guess. If more than one, then moment of whole island is one of the
        guesses if mom1 is within n pixels of one of the maxima. Else dont take
        whole island moment. Instead, find minima on lines connecting all maxima
        and use geometric mean of all minima of a peak as the size of that peak.
        """
        from math import sqrt
        from const import fwsig
        import scipy.ndimage as nd
        import functions as func

        im = isl.image-isl.islmean; mask = isl.mask_active; av = img.clipped_mean; thr1= -1e9
        inipeak, iniposn, im1 = func.get_maxima(im, mask, thr1, isl.shape, beam)
        npeak = len(iniposn)
        gaul = []

        av, stdnew, maxv, maxp, minv, minp = func.arrstatmask(im, mask)
        mom = func.momanalmask_gaus(isl.image-isl.islmean, isl.mask_active, 0, 1.0, True)
        if npeak <= 1:
          g = (float(maxv), int(round(mom[1])), int(round(mom[2])), mom[3]/fwsig, \
                                  mom[4]/fwsig, mom[5])
          gaul.append(g)

        if npeak > 1:                   # markers start from 1=background, watershed starts from 1=background
          watershed, markers = func.watershed(im, mask=isl.mask_active)
          nshed = N.max(markers)-1      # excluding background
          xm, ym = N.transpose([N.where(markers==i) for i in range(1,nshed+2)])[0]
          coords = [c for c in N.transpose([xm,ym])[1:]]
          alldists = [func.dist_2pt(c1, c2) for c1 in coords for c2 in coords if N.any(c1!=c2)]  # has double
          meandist = N.mean(alldists)    # mean dist between all pairs of markers
          compact = []; invmask = []
          for ished in range(nshed):
            shedmask = N.where(watershed==ished+2, False, True) + isl.mask_active # good unmasked pixels = 0
            imm = nd.binary_dilation(~shedmask, N.ones((3,3), int))
            xbad, ybad = N.where((imm==1)*(im>im[xm[ished+1], ym[ished+1]]))
            imm[xbad, ybad] = 0
            invmask.append(imm); x, y = N.where(imm); xcen, ycen = N.mean(x), N.mean(y) # good pixels are now = 1
            dist = func.dist_2pt([xcen, ycen], [xm[ished+1], ym[ished+1]])
            if dist < max(3.0, meandist/4.0):
              compact.append(True)  # if not compact, break source + diffuse
            else:
              compact.append(False)
          if not N.all(compact):
           avsize = []
           ind = N.where(compact)[0]
           for i in ind: avsize.append(N.sum(invmask[i]))
           avsize = sqrt(N.mean(N.array(avsize)))
           for i in range(len(compact)):
             if not compact[i]:                         # make them all compact
               newmask = N.zeros(imm.shape, bool)
               newmask[max(0,xm[i+1]-avsize/2):min(im.shape[0],xm[i+1]+avsize/2), \
                       max(0,ym[i+1]-avsize/2):min(im.shape[1],ym[i+1]+avsize/2)] = True
               invmask[i] = invmask[i]*newmask
          resid = N.zeros(im.shape)                    # approx fit all compact ones
          for i in range(nshed):
            mask1 = ~invmask[i]
            size = sqrt(N.sum(invmask))/fwsig
            xf, yf = coords[i][0], coords[i][1]
            p_ini = [im[xf, yf], xf, yf, size, size, 0.0]
            x, y = N.indices(im.shape)
            p, success = func.fit_gaus2d(im*invmask[i], p_ini, x, y)
            resid = resid + func.gaus_2d(p, x, y)
            gaul.append(p)
          resid = im - resid
          if not N.all(compact):                        # just add one gaussian to fit whole unmasked island
            maxv = N.max(resid)                         # assuming resid has only diffuse emission. can be false
            x, y = N.where(~isl.mask_active); xcen = N.mean(x); ycen = N.mean(y)
            invm = ~isl.mask_active
            #bound = invm - nd.grey_erosion(invm, footprint = N.ones((3,3), int)) # better to use bound for ellipse fitting
            mom = func.momanalmask_gaus(invm, N.zeros(invm.shape, int), 0, 1.0, True)
            g = (maxv, xcen, ycen, mom[3]/fwsig, mom[4]/fwsig, mom[5]-90.)
            gaul.append(g)
            coords.append([xcen, ycen])

        return gaul


    def fit_iter(self, gaul, ng1, fcn, dof, beam, thr, iter, inifit, ngmax, verbose=1):
        """One round of fitting

        Parameters:
        gaul : list of initial gaussians
        fcn  : MGFunction object
        dof  : maximal number of fitted parameters
        beam : initial shape for newly added gaussians
               [bmaj, bmin, bpa] in pixels
        thr  : peak threshold for adding more gaussians
        verbose: whether to print fitting progress information
        """
        from _cbdsm import lmder_fit, dn2g_fit, dnsg_fit
#         global bar
        fit = lmder_fit
        beam = list(beam)

        ### first drop-in initial gaussians
        ### no error-checking here, they MUST fit
        fcn.reset()
        for ig in range(ng1):
          g = gaul[ig]
          self.add_gaussian(fcn, g, dof)

        ### do a round of fitting if any initials were provided
        fitok = True
        if len(gaul) != 0:
          fitok = fit(fcn, final=0, verbose=verbose)

        ### iteratively add gaussians while there are high peaks
        ### in the image and fitting converges
        while fitok:
#           if bar.started: bar.spin()
          peak, coords = fcn.find_peak()
          if peak < thr:  ### no good peaks left
              break
          if len(fcn.parameters) < ngmax and iter == 1 and inifit == 'default' and len(gaul) >= ng1+1:
             ng1 = ng1 + 1
             g = gaul[ng1-1]
          else:
            if len(fcn.parameters) < ngmax and inifit in ['simple', 'nobeam']:
              g = [peak, coords[0], coords[1]] + beam
            else:
              break
          fitok &= self.add_gaussian(fcn, g, dof)

          fitok &= fit(fcn, final=0, verbose=verbose)

        ### and one last fit with higher precision
        ### make sure we return False when fitok==False due to lack
        ### of free parameters
        fitok &= fit(fcn, final=1, verbose=verbose)
        return fitok

    def add_gaussian(self, fcn, parameters, dof):
        """Try adding one more gaussian to fcn object.
        It's trying to reduce number of fitted parameters if
        there is not enough DoF left.

        Parameters:
        fcn: MGFunction object
        parameters: initial values for gaussian parameters
        dof: total possible number of fitted parameters
        """
        from _cbdsm import Gtype

        gtype = (Gtype.g3 if fcn.fitted_parameters() + 3 <= dof else None)
        gtype = (Gtype.g6 if fcn.fitted_parameters() + 6 <= dof else gtype)

        if gtype:
            fcn.add_gaussian(gtype, parameters)
            return True
        else:
            return False

    def flag_gaussians(self, gaul, opts, beam, thr, peak, shape, isl_mask, isl_image, size):
        """Flag gaussians according to some rules.
        Splits list of gaussian parameters in 2, where the first
        one is a list of parameters for accepted gaussians, and
        the second one is a list of pairs (flag, parameters) for
        flagged gaussians.

        Parameters:
        gaul: input list of gaussians
        opts: Opts object to extract flagging parameters from
        beam: beam shape
        thr: threshold for pixels with signal
        peak: peak data value in the current island
        shape: shape of the current island
        isl_mask: island mask
        """
        good = []
        bad  = []
        for g in gaul:

            flag = self._flag_gaussian(g, beam, thr, peak, shape, opts, isl_mask, isl_image, size)
            if flag:
                bad.append((flag, g))
            else:
                good.append(g)

        return good, bad

    def _flag_gaussian(self, g, beam, thr, peak, shape, opts, mask, image, size_bms):
        """The actual flagging routine. See above for description.
        """
        from math import sqrt, sin, cos, log, pi
        from const import fwsig
        import functions as func
        import scipy.ndimage as nd

        A, x1, x2, s1, s2, th = g
        s1, s2 = map(abs, [s1, s2])
        flag = 0
        if N.any(N.isnan(g)) or s1 == 0.0 or s2 == 0.0:
            return -1

        if s1 < s2:   # s1 etc are sigma
          ss1=s2; ss2=s1; th1 = divmod(th+90.0, 180)[1]
        else:
          ss1=s1; ss2=s2; th1 = divmod(th, 180)[1]
        th1 = th1/180.0*pi
        if ss1 > 1e4 and ss2 > 1e4:
          xbox = 1e9; ybox = 1e9
        else:
          xbox = 2.0*(abs(ss1*cos(th1)*cos(th1))+abs(ss2*ss2/ss1*sin(th1)*sin(th1)))/ \
                 sqrt(cos(th1)*cos(th1)+ss2*ss2/ss1/ss1*sin(th1)*sin(th1))
          ybox = 2.0*(abs(ss1*sin(th1)*sin(th1))+abs(ss2*ss2/ss1*cos(th1)*cos(th1)))/ \
                 sqrt(sin(th1)*sin(th1)+ss2*ss2/ss1/ss1*cos(th1)*cos(th1))

        ### now check all conditions
        border = opts.flag_bordersize
        x1ok = True
        x2ok = True
        flagmax = False
        if A < opts.flag_minsnr*thr: flag += 1
        if A > opts.flag_maxsnr*peak:
            flag += 2
            flagmax = True
        if x1 - border < 0 or x1 + border + 1 > shape[0]:
            flag += 4
            x1ok = False
        if x2 - border < 0 or x2 + border + 1 > shape[1]:
            flag += 8
            x2ok = False
        if x1ok and x2ok:
            if not flagmax:
                # Check image value at Gaussian center
                im_val_at_cen = nd.map_coordinates(image, [N.array([x1]), N.array([x2])])
                if A > opts.flag_maxsnr*im_val_at_cen:
                   flag += 2
            borx1_1 = x1 - border
            if borx1_1 < 0: borx1_1 = 0
            borx1_2 = x1 + border + 1
            if borx1_2 > shape[0]: borx1_2 = shape[0]
            if N.any(mask[borx1_1:borx1_2, x2]):
                flag += 4
            borx2_1 = x2 - border
            if borx2_1 < 0: borx2_1 = 0
            borx2_2 = x2 + border + 1
            if borx2_2 > shape[1]: borx2_2 = shape[1]
            if N.any(mask[x1, borx2_1:borx2_2]):
                flag += 8
        if xbox > opts.flag_maxsize_isl*shape[0]: flag += 16
        if ybox > opts.flag_maxsize_isl*shape[1]: flag += 32
        if s1*s2 > opts.flag_maxsize_bm*beam[0]*beam[1]: flag += 64
        if opts.flag_smallsrc:
          if s1*s2 < opts.flag_minsize_bm*beam[0]*beam[1]: flag += 128
        if not opts.flag_smallsrc:
                if s1*s2 == 0.: flag += 128

        if size_bms > 30.0:
            # Only check if island is big enough, as this flagging step
            # is unreliable for small islands. size_bms is size of island
            # in number of beam areas
            ellx, elly = func.drawellipse([A, x1, x2, s1*opts.flag_maxsize_fwhm,
                                           s2*opts.flag_maxsize_fwhm, th])
            pt1 = [N.min(ellx), elly[N.argmin(ellx)]]
            pt2 = [ellx[N.argmax(elly)], N.max(elly)]
            pt3 = [N.max(ellx), elly[N.argmax(ellx)]]
            pt4 = [ellx[N.argmin(elly)], N.min(elly)]
            extremes = [pt1, pt2, pt3, pt4]
            for pt in extremes:
                if pt[0] < 0 or pt[0] >= shape[0] or pt[1] < 0 or pt[1] >= shape[1]:
                    flag += 256
                    break
                elif mask[tuple(pt)]:
                    flag += 256
                    break

        return flag

    def fixup_gaussian(self, isl, gaussian):
        """Normalize parameters by adjusting them to the
        proper image coordinates and ensuring that all of
        the implicit conventions (such as bmaj >= bmin) are met.
        """
        np = list(gaussian)

        ### update to the image coordinates
        np[1] += isl.origin[0]
        np[2] += isl.origin[1]

        ### shape values should be positive
        np[3] = abs(np[3])
        np[4] = abs(np[4])

        ### first extent is major
        if np[3] < np[4]:
            np[3:5] = np[4:2:-1]
            np[5] += 90

        ### clip position angle
        np[5] = divmod(np[5], 180)[1]

        return np

    def check_for_1d_gaussians(self, img):
        """Check for Gaussians with deconvolved sizes of 0 for one axis only."""
        n1d = 0
        ng = 0
        for g in img.gaussians:
            ng += 1
            dsize = g.deconv_size_sky
            if (dsize[0] == 0 and dsize[1] > 0) or (dsize[0] > 0 and dsize[1] == 0):
                n1d += 1
        if ng > 0:
            return float(n1d)/float(ng)
        else:
            return 0.0

def find_bbox(thresh, g):
    """Calculate bounding box for gaussian.

    This function calculates size of the box for evaluating
    gaussian, so that value of gaussian is smaller than threshold
    outside of the box.

    Parameters:
    thres: threshold
    g: Gaussian object or list of paramters
    """

    from math import ceil, sqrt, log
    if isinstance(g, list):
        A = g[0]
        S = g[3]
    else:
        A = g.peak_flux
        S = g.size_pix[0]
    if A == 0.0:
        return ceil(S*1.5)
    if thresh/A >= 1.0 or thresh/A <= 0.0:
        return ceil(S*1.5)
    return ceil(S*sqrt(-2*log(thresh/A)))


from image import *

class Gaussian(object):
    """Instances of this class are used to store information about
    extracted gaussians in a structured way.
    """
    gaussian_idx = Int(doc="Serial number of the gaussian within island")
    gaus_num = Int(doc="Serial number of the gaussian for the image", colname='Gaus_id')
    island_id   = Int(doc="Serial number of the island", colname='Isl_id')
    flag        = Int(doc="Flag associated with gaussian", colname='Flag')
    parameters  = List(Float(), doc="Raw gaussian parameters")
    total_flux  = Float(doc="Total flux density, Jy", colname='Total_flux', units='Jy')
    total_fluxE = Float(doc="Total flux density error, Jy", colname='E_Total_flux',
                        units='Jy')
    peak_flux   = Float(doc="Peak flux density/beam, Jy/beam", colname='Peak_flux',
                        units='Jy/beam')
    peak_fluxE  = Float(doc="Peak flux density/beam error, Jy/beam",
                        colname='E_Peak_flux', units='Jy/beam')
    centre_sky  = List(Float(), doc="Sky coordinates of gaussian centre",
                       colname=['RA', 'DEC'], units=['deg', 'deg'])
    centre_skyE = List(Float(), doc="Error on sky coordinates of gaussian centre",
                       colname=['E_RA', 'E_DEC'], units=['deg', 'deg'])
    centre_pix  = List(Float(), doc="Pixel coordinates of gaussian centre",
                       colname=['Xposn', 'Yposn'], units=['pix', 'pix'])
    centre_pixE = List(Float(), doc="Error on pixel coordinates of gaussian centre",
                       colname=['E_Xposn', 'E_Yposn'], units=['pix', 'pix'])
    size_sky   = List(Float(), doc="Shape of the gaussian FWHM, PA, deg",
                      colname=['Maj', 'Min', 'PA'], units=['deg', 'deg',
                      'deg'])
    size_skyE  = List(Float(), doc="Error on shape of the gaussian FWHM, PA, deg",
                      colname=['E_Maj', 'E_Min', 'E_PA'], units=['deg', 'deg',
                      'deg'])
    deconv_size_sky = List(Float(), doc="Deconvolved shape of the gaussian FWHM, PA, deg",
                      colname=['DC_Maj', 'DC_Min', 'DC_PA'], units=['deg', 'deg',
                      'deg'])
    deconv_size_skyE = List(Float(), doc="Error on deconvolved shape of the gaussian FWHM, PA, deg",
                      colname=['E_DC_Maj', 'E_DC_Min', 'E_DC_PA'], units=['deg', 'deg',
                      'deg'])
    size_pix   = List(Float(), doc="Shape of the gaussian FWHM, pixel units")
    size_pixE  = List(Float(), doc="Error on shape of the gaussian FWHM, pixel units")
    rms        = Float(doc="Island rms Jy/beam", colname='Isl_rms', units='Jy/beam')
    mean       = Float(doc="Island mean Jy/beam", colname='Isl_mean', units='Jy/beam')
    total_flux_isl = Float(doc="Island total flux from sum of pixels", colname='Isl_Total_flux', units='Jy')
    total_flux_islE = Float(doc="Error on island total flux from sum of pixels", colname='E_Isl_Total_flux', units='Jy')
    gresid_rms = Float(doc="Island rms in Gaussian residual image", colname='Resid_Isl_rms', units='Jy/beam')
    gresid_mean= Float(doc="Island mean in Gaussian residual image", colname='Resid_Isl_mean', units='Jy/beam')
    sresid_rms = Float(doc="Island rms in Shapelet residual image", colname='Resid_Isl_rms', units='Jy/beam')
    sresid_mean= Float(doc="Island mean in Shapelet residual image", colname='Resid_Isl_mean', units='Jy/beam')
    jlevel     = Int(doc="Wavelet number to which Gaussian belongs", colname='Wave_id')

    def __init__(self, img, gaussian, isl_idx, g_idx, flag=0):
        """Initialize Gaussian object from fitting data

        Parameters:
        img: PyBDSM image object
        gaussian: 6-tuple of fitted numbers
        isl_idx: island serial number
        g_idx: gaussian serial number
        flag: flagging (if any)
        """
        import functions as func
        from const import fwsig
        import numpy as N

        self.gaussian_idx = g_idx
        self.gaus_num = 0 # stored later
        self.island_id = isl_idx
        self.jlevel = img.j
        self.flag = flag
        self.parameters = gaussian

        p = gaussian
        self.peak_flux = p[0]
        self.centre_pix = p[1:3]
        size = p[3:6]
        if func.approx_equal(size[0], img.pixel_beam[0]*1.1) and \
                func.approx_equal(size[1], img.pixel_beam[1]) and \
                func.approx_equal(size[2], img.pixel_beam[2]):
            # Check whether fitted Gaussian is just the distorted pixel beam
            # given as an initial guess. If so, reset the size to the
            # undistorted beam.
            size = img.pixel_beam
        size = func.corrected_size(size)  # gives fwhm and P.A.
        self.size_pix = size # FWHM in pixels and P.A. CCW from +y axis
        self.size_sky = img.pix2beam(size, self.centre_pix) # FWHM in degrees and P.A. CCW from North

        # Check if this is a wavelet image. If so, use orig_pixel_beam
        # for flux calculation, as pixel_beam has been altered to match
        # the wavelet scale.
        if img.waveletimage:
            pixel_beam = img.orig_pixel_beam
        else:
            pixel_beam = img.pixel_beam
        bm_pix = N.array([pixel_beam[0]*fwsig, pixel_beam[1]*fwsig, pixel_beam[2]])
        tot = p[0]*size[0]*size[1]/(bm_pix[0]*bm_pix[1])

        if flag == 0:
          errors = func.get_errors(img, p+[tot], img.islands[isl_idx].rms)
          self.centre_sky = img.pix2sky(p[1:3])
        else:
          errors = [0]*7
          self.centre_sky = [0., 0.]
        self.total_flux = tot
        self.total_fluxE = errors[6]

        self.peak_fluxE = errors[0]
        self.total_fluxE = errors[6]
        self.centre_pixE = errors[1:3]
        self.centre_skyE = img.pix2coord(errors[1:3])
        self.size_pixE = errors[3:6]
        self.size_skyE = img.pix2beam(errors[3:6], self.centre_pix)
        self.rms = img.islands[isl_idx].rms
        self.mean = img.islands[isl_idx].mean
        self.total_flux_isl = img.islands[isl_idx].total_flux
        self.total_flux_islE = img.islands[isl_idx].total_fluxE

        # func.deconv, based on AIPS DECONV.FOR, gives a lot of
        # 1-D Gaussians. The Miriad algorithm in func.deconv2 does
        # a much better job in this regard, so use it instead. Note
        # that for resolved sources, the two algorithms give the
        # same answer. For now, errors on the deconvolved parameters
        # are just set to those of the undeconvolved ones.
        if flag == 0:
            gaus_dc, err = func.deconv2(bm_pix, size)
            self.deconv_size_sky = img.pix2beam(gaus_dc, self.centre_pix)
            self.deconv_size_skyE  = img.pix2beam(errors[3:6], self.centre_pix)
        else:
            self.deconv_size_sky = [0., 0., 0.]
            self.deconv_size_skyE  = [0., 0., 0.]


### Insert attributes into Island class
from islands import Island
Island.gaul = List(tInstance(Gaussian), doc="List of extracted gaussians")
Island.fgaul= List(tInstance(Gaussian),
                   doc="List of extracted (flagged) gaussians")
