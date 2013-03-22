
import numpy as N
from image import *
import mylogger
from copy import deepcopy as cp
from . import has_pl
if has_pl:
    import matplotlib.pyplot as pl
import scipy
import scipy.signal as S
import _cbdsm
import functions as func
import _pytesselate as _pytess
import shapelets as sh
from scipy.optimize import leastsq
import nat
from math import *
import statusbar
from const import fwsig
import multi_proc as mp
import itertools


class Op_psf_vary(Op):
    """Computes variation of psf across the image """

    def __call__(self, img):

      if img.opts.psf_vary_do:
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Psf_Vary")
        mylogger.userinfo(mylog, '\nEstimating PSF variations')
        opts = img.opts
        dir = img.basedir + '/misc/'
        plot = False # debug figures
        image = img.ch0

        over = 2
        generators = opts.psf_generators; nsig = opts.psf_nsig; kappa2 = opts.psf_kappa2
        snrtop = opts.psf_snrtop; snrbot = opts.psf_snrbot; snrcutstack = opts.psf_snrcutstack
        gencode = opts.psf_gencode; primarygen = opts.psf_primarygen; itess_method = opts.psf_itess_method
        tess_sc = opts.psf_tess_sc; tess_fuzzy= opts.psf_tess_fuzzy
        bright_snr_cut = opts.psf_high_snr
        s_only = opts.psf_stype_only
        if opts.psf_snrcut < 5.0:
            mylogger.userinfo(mylog, "Value of psf_snrcut too low; increasing to 5")
            snrcut = 5.0
        else:
            snrcut = opts.psf_snrcut
        img.psf_snrcut = snrcut
        if opts.psf_high_snr != None:
            if opts.psf_high_snr < 10.0:
                mylogger.userinfo(mylog, "Value of psf_high_snr too low; increasing to 10")
                high_snrcut = 10.0
            else:
                high_snrcut = opts.psf_high_snr
        else:
            high_snrcut = opts.psf_high_snr
        img.psf_high_snr = high_snrcut

        wtfns=['unity', 'roundness', 'log10', 'sqrtlog10']
        if 0 <= itess_method < 4: tess_method=wtfns[itess_method]
        else: tess_method='unity'

        ### now put all relevant gaussian parameters into a list
        ngaus = img.ngaus
        nsrc = img.nsrc
        num = N.zeros(nsrc, int)
        peak = N.zeros(nsrc)
        xc = N.zeros(nsrc)
        yc = N.zeros(nsrc)
        bmaj = N.zeros(nsrc)
        bmin = N.zeros(nsrc)
        bpa = N.zeros(nsrc)
        code = N.array(['']*nsrc);
        rms = N.zeros(nsrc)
        src_id_list = []
        for i, src in enumerate(img.sources):
            src_max = 0.0
            for gmax in src.gaussians:
                # Take only brightest Gaussian per source
                if gmax.peak_flux > src_max:
                    src_max = gmax.peak_flux
                    g = gmax
            num[i] = i
            peak[i] = g.peak_flux
            xc[i] = g.centre_pix[0]
            yc[i] = g.centre_pix[1]
            bmaj[i] = g.size_pix[0]
            bmin[i] = g.size_pix[1]
            bpa[i] = g.size_pix[2]
            code[i] = img.sources[g.source_id].code
            rms[i] = img.islands[g.island_id].rms
        gauls = (num, peak, xc, yc, bmaj, bmin, bpa, code, rms)
        tr_gauls = self.trans_gaul(gauls)

        # takes gaussians with code=S and snr > snrcut.
        if s_only:
            tr = [n for n in tr_gauls if n[1]/n[8]>snrcut and n[7] == 'S']
        else:
            tr = [n for n in tr_gauls if n[1]/n[8]>snrcut]
        g_gauls = self.trans_gaul(tr)

        # computes statistics of fitted sizes. Same as psfvary_fullstat.f in fBDSM.
        bmaj_a, bmaj_r, bmaj_ca, bmaj_cr, ni = _cbdsm.bstat(bmaj, None, nsig)
        bmin_a, bmin_r, bmin_ca, bmin_cr, ni = _cbdsm.bstat(bmin, None, nsig)
        bpa_a, bpa_r, bpa_ca, bpa_cr, ni = _cbdsm.bstat(bpa, None, nsig)

        # get subset of sources deemed to be unresolved. Same as size_ksclip_wenss.f in fBDSM.
        flag_unresolved = self.get_unresolved(g_gauls, img.beam, nsig, kappa2, over, img.psf_high_snr, plot)

        # see how much the SNR-weighted sizes of unresolved sources differ from the synthesized beam.
        wtsize_beam_snr = self.av_psf(g_gauls, img.beam, flag_unresolved)

        # filter out resolved sources
        tr_gaul = self.trans_gaul(g_gauls)
        tr = [n for i, n in enumerate(tr_gaul) if flag_unresolved[i]]
        g_gauls = self.trans_gaul(tr)
        mylogger.userinfo(mylog, 'Number of unresolved sources', str(len(g_gauls[0])))

        # get a list of voronoi generators. vorogenS has values (and not None) if generators='field'.
        vorogenP, vorogenS = self.get_voronoi_generators(g_gauls, generators, gencode, snrcut, snrtop, snrbot, snrcutstack)
        mylogger.userinfo(mylog, 'Number of generators for PSF variation', str(len(vorogenP[0])))
        if len(vorogenP[0]) < 3:
            mylog.warning('Insufficient number of generators')
            return

        mylogger.userinfo(mylog, 'Tesselating image')
        # group generators into tiles
        tile_prop = self.edit_vorogenlist(vorogenP, frac=0.9)

        # tesselate the image
        volrank, vorowts = self.tesselate(vorogenP, vorogenS, tile_prop, tess_method, tess_sc, tess_fuzzy, \
                  generators, gencode, image.shape)
        if opts.output_all:
            func.write_image_to_file(img.use_io, img.imagename + '.volrank.fits', volrank, img, dir)

        tile_list, tile_coord, tile_snr = tile_prop
        ntile = len(tile_list)
        bar = statusbar.StatusBar('Determining PSF variation ............... : ', 0, ntile)
        mylogger.userinfo(mylog, 'Number of tiles for PSF variation', str(ntile))

        # For each tile, calculate the weighted averaged psf image. Also for all the sources in the image.
        cdelt = list(img.wcs_obj.acdelt[0:2])
        factor=3.
        psfimages, psfcoords, totpsfimage, psfratio, psfratio_aper = self.psf_in_tile(image, img.beam, g_gauls, \
                   cdelt, factor, snrcutstack, volrank, tile_prop, plot, img)
        npsf = len(psfimages)

        if opts.psf_use_shap:
            # use totpsfimage to get beta, centre and nmax for shapelet decomposition. Use nmax=5 or 6
            mask=N.zeros(totpsfimage.shape, dtype=bool)
            (m1, m2, m3)=func.moment(totpsfimage, mask)
            betainit=sqrt(m3[0]*m3[1])*2.0  * 1.4
            tshape = totpsfimage.shape
            cen = N.array(N.unravel_index(N.argmax(totpsfimage), tshape))+[1,1]
            cen = tuple(cen)
            nmax = 12
            basis = 'cartesian'
            betarange = [0.5,sqrt(betainit*max(tshape))]
            beta, error  = sh.shape_varybeta(totpsfimage, mask, basis, betainit, cen, nmax, betarange, plot)
            if error == 1: print '  Unable to find minimum in beta'

            # decompose all the psf images using the beta from above
            nmax=12; psf_cf=[]
            for i in range(npsf):
                psfim = psfimages[i]
                cf = sh.decompose_shapelets(psfim, mask, basis, beta, cen, nmax, mode='')
                psf_cf.append(cf)
                if img.opts.quiet == False:
                    bar.increment()
            bar.stop()

            # transpose the psf image list
            xt, yt = N.transpose(tile_coord)
            tr_psf_cf = N.transpose(N.array(psf_cf))

            # interpolate the coefficients across the image. Ok, interpolate in scipy for
            # irregular grids is crap. doesnt even pass through some of the points.
            # for now, fit polynomial.
            compress = 100.0
            x, y = N.transpose(psfcoords)
            if len(x) < 3:
                mylog.warning('Insufficient number of tiles to do interpolation of PSF variation')
                return

            psf_coeff_interp, xgrid, ygrid = self.interp_shapcoefs(nmax, tr_psf_cf, psfcoords, image.shape, \
                     compress, plot)

            psfshape = psfimages[0].shape
            skip = 5
            aa = self.create_psf_grid(psf_coeff_interp, image.shape, xgrid, ygrid, skip, nmax, psfshape, \
                 basis, beta, cen, totpsfimage, plot)
            img.psf_images = aa
        else:
            if ntile < 4:
                mylog.warning('Insufficient number of tiles to do interpolation of PSF variation')
                return
            else:
                # Fit stacked PSFs with Gaussians and measure aperture fluxes
                bm_pix = N.array([img.pixel_beam()[0]*fwsig, img.pixel_beam()[1]*fwsig, img.pixel_beam()[2]])
                psf_maj = N.zeros(npsf)
                psf_min = N.zeros(npsf)
                psf_pa = N.zeros(npsf)
                if img.opts.quiet == False:
                    bar.start()
                for i in range(ntile):
                    psfim = psfimages[i]
                    mask = N.zeros(psfim.shape, dtype=bool)
                    x_ax, y_ax = N.indices(psfim.shape)
                    maxv = N.max(psfim)
                    p_ini = [maxv, (psfim.shape[0]-1)/2.0*1.1, (psfim.shape[1]-1)/2.0*1.1, bm_pix[0]/fwsig*1.3,
                             bm_pix[1]/fwsig*1.1, bm_pix[2]*2]
                    para, ierr = func.fit_gaus2d(psfim, p_ini, x_ax, y_ax, mask)
                    ### first extent is major
                    if para[3] < para[4]:
                        para[3:5] = para[4:2:-1]
                        para[5] += 90
                    ### clip position angle
                    para[5] = divmod(para[5], 180)[1]

                    psf_maj[i] = para[3]
                    psf_min[i] = para[4]
                    posang = para[5]
                    while posang >= 180.0:
                        posang -= 180.0
                    psf_pa[i] = posang

                    if img.opts.quiet == False:
                        bar.increment()
                bar.stop()

                # Interpolate Gaussian parameters
                if img.aperture == None:
                    psf_maps = [psf_maj, psf_min, psf_pa, psfratio]
                else:
                    psf_maps = [psf_maj, psf_min, psf_pa, psfratio, psfratio_aper]
                nimgs = len(psf_maps)
                bar = statusbar.StatusBar('Interpolating PSF images ................ : ', 0, nimgs)
                if img.opts.quiet == False:
                    bar.start()
                map_list = mp.parallel_map(func.eval_func_tuple,
                    itertools.izip(itertools.repeat(self.interp_prop),
                    psf_maps, itertools.repeat(psfcoords),
                    itertools.repeat(image.shape)), numcores=opts.ncores,
                    bar=bar)
                if img.aperture == None:
                    psf_maj_int, psf_min_int, psf_pa_int, psf_ratio_int = map_list
                else:
                    psf_maj_int, psf_min_int, psf_pa_int, psf_ratio_int, psf_ratio_aper_int = map_list

                # Smooth if desired
                if img.opts.psf_smooth != None:
                    sm_scale = img.opts.psf_smooth / img.pix2beam([1.0, 1.0, 0.0])[0] / 3600.0 # pixels
                    if img.opts.aperture == None:
                        psf_maps = [psf_maj_int, psf_min_int, psf_pa_int, psf_ratio_int]
                    else:
                        psf_maps = [psf_maj_int, psf_min_int, psf_pa_int, psf_ratio_int, psf_ratio_aper_int]
                    nimgs = len(psf_maps)
                    bar = statusbar.StatusBar('Smoothing PSF images .................... : ', 0, nimgs)
                    if img.opts.quiet == False:
                        bar.start()
                    map_list = mp.parallel_map(func.eval_func_tuple,
                        itertools.izip(itertools.repeat(self.blur_image),
                        psf_maps, itertools.repeat(sm_scale)), numcores=opts.ncores,
                        bar=bar)
                    if img.aperture == None:
                        psf_maj_int, psf_min_int, psf_pa_int, psf_ratio_int = map_list
                    else:
                        psf_maj_int, psf_min_int, psf_pa_int, psf_ratio_int, psf_ratio_aper_int = map_list

                # Make sure all smoothed, interpolated images are ndarrays
                psf_maj_int = N.array(psf_maj_int)
                psf_min_int = N.array(psf_min_int)
                psf_pa_int = N.array(psf_pa_int)
                psf_ratio_int = N.array(psf_ratio_int)
                if img.aperture == None:
                    psf_ratio_aper_int = N.zeros(psf_maj_int.shape)
                else:
                    psf_ratio_aper_int = N.array(psf_ratio_aper_int)

                # Store interpolated images. The major and minor axis images are
                # the sigma in units of arcsec, the PA image in units of degrees east of
                # north, the ratio images in units of 1/beam.
                img.psf_vary_maj = psf_maj_int * img.pix2beam([1.0, 1.0, 0.0])[0] * 3600.0 # sigma in arcsec
                img.psf_vary_min = psf_min_int * img.pix2beam([1.0, 1.0, 0.0])[0] * 3600.0 # sigma in arcsec
                img.psf_vary_pa = psf_pa_int
                img.psf_vary_ratio = psf_ratio_int # in 1/beam
                img.psf_vary_ratio_aper = psf_ratio_aper_int # in 1/beam

                # Blank with NaNs if needed
                mask = img.mask
                if isinstance(mask, N.ndarray):
                    pix_masked = N.where(mask == True)
                    img.psf_vary_maj[pix_masked] = N.nan
                    img.psf_vary_min[pix_masked] = N.nan
                    img.psf_vary_pa[pix_masked] = N.nan
                    img.psf_vary_ratio[pix_masked] = N.nan
                    img.psf_vary_ratio_aper[pix_masked] = N.nan

                if opts.output_all:
                    func.write_image_to_file(img.use_io, img.imagename + '.psf_vary_maj.fits', img.psf_vary_maj*fwsig, img, dir)
                    func.write_image_to_file(img.use_io, img.imagename + '.psf_vary_min.fits', img.psf_vary_min*fwsig, img, dir)
                    func.write_image_to_file(img.use_io, img.imagename + '.psf_vary_pa.fits', img.psf_vary_pa, img, dir)
                    func.write_image_to_file(img.use_io, img.imagename + '.psf_vary_ratio.fits', img.psf_vary_ratio, img, dir)
                    func.write_image_to_file(img.use_io, img.imagename + '.psf_vary_ratio_aper.fits', img.psf_vary_ratio_aper, img, dir)

                # Loop through source and Gaussian lists and deconvolve the sizes using appropriate beam
                bar2 = statusbar.StatusBar('Correcting deconvolved source sizes ..... : ', 0, img.nsrc)
                if img.opts.quiet == False:
                    bar2.start()
                for src in img.sources:
                    src_pos = img.sky2pix(src.posn_sky_centroid)
                    src_pos_int = (int(src_pos[0]), int(src_pos[1]))
                    gaus_c = img.gaus2pix(src.size_sky, src.posn_sky_centroid)
                    gaus_bm = [psf_maj_int[src_pos_int]*fwsig, psf_min_int[src_pos_int]*fwsig, psf_pa_int[src_pos_int]*fwsig]
                    gaus_dc, err = func.deconv2(gaus_bm, gaus_c)
                    src.deconv_size_sky = img.pix2gaus(gaus_dc, src_pos)
                    src.deconv_size_skyE = [0.0, 0.0, 0.0]
                    for g in src.gaussians:
                        gaus_c = img.gaus2pix(g.size_sky, src.posn_sky_centroid)
                        gaus_dc, err = func.deconv2(gaus_bm, gaus_c)
                        g.deconv_size_sky = img.pix2gaus(gaus_dc, g.centre_pix)
                        g.deconv_size_skyE = [0.0, 0.0, 0.0]
                        if img.opts.quiet == False:
                            bar2.spin()
                    if img.opts.quiet == False:
                        bar2.increment()
                bar2.stop()
        img.completed_Ops.append('psf_vary')

##################################################################################################

    def trans_gaul(self, q):
        " transposes a tuple of .gaul values "
        y=[]
        for i in range(len(q[0])):
            elem=[]
            for j in range(len(q)):
                elem.append(q[j][i])
            y.append(elem)
        return y

##################################################################################################

    def bindata(self, over,num): #ptpbin,nbin,ptplastbin, same as get_bins in fBDSM.

        if num <= 100: ptpbin=num/5
        if num > 100: ptpbin=num/10
        if num > 1000: ptpbin=num/20
        if ptpbin % 2 == 1: ptpbin=ptpbin+1
        if num < 10: ptpbin=num
        nbin=(num-ptpbin)/(ptpbin/over)+1
        ptplastbin=(num-1)-(nbin-1)*ptpbin/over
        nbin=nbin+1

        return ptpbin, nbin, ptplastbin

##################################################################################################
    def bin_and_stats_ny(self, x,y,over,ptpbin,nbin,ptplastbin,nsig):
        import math

        n1=N.array(range(nbin))+1   # bin number
        n2=N.array([ptpbin]*nbin); n2[nbin-2]=ptplastbin; n2[nbin-1]=ptpbin/over
        n3=N.array([ptpbin]*nbin, dtype=float); n3[nbin-1]=float(over)*(len(x)-ptpbin/2)/(nbin-1)
        xval=N.zeros(nbin)
        meany=N.zeros(nbin); stdy=N.zeros(nbin); mediany=N.zeros(nbin)
        for i in range(nbin):
            lb=round(1+(n1[i]-1)*n3[i]/over+(1-1))-1 # -1 for python indexing
            ub=round(1+(n1[i]-1)*n3[i]/over+(n2[i]-1))-1 # -1 for python indexing
            x1=x[lb:ub+1]; y1=y[lb:ub+1]

            # do calcmedianclip2vec.f for code=YYN
            nout=100; niter=0
            while nout>0 and niter<6:
              med1=N.median(y1[:])
              med2=10.**(N.median(N.log10(x1[:])))
              medstd=0    # calcmedianstd.f
              for j in y1: medstd += (j-med1)*(j-med1)
              medstd=math.sqrt(medstd/len(y1))        #
              av1=N.mean(y1); std1=func.std(y1)
              av2=N.mean(x1); std2=func.std(x1)
              # get_medianclip_vec2
              z=N.transpose([x1, y1])
              z1=N.transpose([n for n in z if abs(n[1]-med1)<=nsig*medstd])
              nout=len(x1)-len(z1[0])
              x1=z1[0]; y1=z1[1];
              niter+=1
            xval[i]=med2;
            meany[i]=av1; stdy[i]=std1; mediany[i]=med1

        if stdy[nbin-1]/mediany[nbin-1] > stdy[nbin-2]/mediany[nbin-2]:
           stdy[nbin-1]=stdy[nbin-2]/mediany[nbin-2]*mediany[nbin-1]
        return xval, meany, stdy, mediany

##################################################################################################
    def LM_fit(self, x, y, err, funct, order=0):
        if funct == func.poly:
           p0=N.array([y[N.argmax(x)]] + [0]*order)
        if funct == func.wenss_fit:
           p0=N.array([y[N.argmax(x)]] + [1.])
        res=lambda p, x, y, err: (y-funct(p, x))/err
        (p, flag)=leastsq(res, p0, args=(x, y, err))
        return p

##################################################################################################

    def fit_bins_func(self, x,y,over,ptpbin,nbin,ptplastbin,nsig):  # sub_size_ksclip
        import math

        (xval,meany,stdy,medy)=self.bin_and_stats_ny(x,y,over,ptpbin,nbin,ptplastbin,nsig)
        yfit=stdy/medy
        err=N.array([1.]*nbin)
        err[nbin-2]=err[0]*math.sqrt(1.0*ptpbin/ptplastbin)
        err[nbin-1]=err[0]*math.sqrt(1.0*ptpbin*over/ptplastbin)

        i=0
        while i<nbin-4 and (N.all(N.sort(yfit[i:i+4])[::-1] == yfit[i:i+4]) == False):
           i+=1
        if i==nbin-4: sind=0
        else: sind=i-1
        if sind < 1:
            sind = 0
        if sind > 0.25*nbin:
            sind=int(round(0.25*nbin))-1

        s_c=self.LM_fit(xval[sind:],yfit[sind:],err[sind:], func.wenss_fit)

        err[:]=1.
        s_cm=self.LM_fit(N.log10(xval),medy,err,func.poly, order=1)
        if len(xval) >= 3:
            s_dm=self.LM_fit(N.log10(xval),medy,err,func.poly, order=2)
        else:
            s_dm = (N.array([s_cm[0], s_cm[1], 0.0]), 0)

        if ptpbin<75: s_dm=N.append(s_cm[:], [0.])
        return s_c, s_dm

##################################################################################################
    def get_unresolved(self, g_gauls, beam, nsig, kappa2, over, bright_snr_cut=20.0, plot=False):
        """"Gets subset of unresolved sources

        Also flags as unresolved all sources with SNRs above
        bright_cut_snr, since fitting below is unreliable for bright
        sources.
        """

        num=len(g_gauls[0])
        b1=N.asarray(g_gauls[4])/(beam[0]*3600.)
        b2=N.asarray(g_gauls[5])/(beam[1]*3600.)
        s1=N.asarray(g_gauls[1])/N.array(g_gauls[8])
        snr=N.array(s1)
        index=snr.argsort()
        snr=snr[index]
        nmaj=N.array(b1)[index]
        nmin=N.array(b2)[index]

#         if plot: pl.figure()
        f_sclip=N.zeros((2,num), dtype=bool)
        for idx, nbeam in enumerate([nmaj, nmin]):
          xarr=N.copy(snr)
          yarr=N.copy(nbeam)
          niter=0; nout=num; noutold=nout*2
          while niter<10 and nout >0.75*num:
            (ptpbin, nbin, ptplastbin)=self.bindata(over,nout)    # get_bins in fBDSM
            (s_c,s_dm) = self.fit_bins_func(xarr,yarr,over,ptpbin,nbin,ptplastbin,nsig)  # size_ksclip_wenss in fBDSM
            noutold = len(xarr)
            z = N.transpose([xarr, yarr, s_dm[0]+s_dm[1]*N.log10(xarr)+s_dm[2]*(N.log10(xarr)**2.),  \
                N.sqrt(s_c[0]*s_c[0]+s_c[1]*s_c[1]/(xarr*xarr)) ])
            z1 = N.transpose([n for n in z if abs(n[1]-n[2])/(n[2]*n[3])<kappa2])  # sub_size_wenss_getnum in fBDSM
            nout = len(z1[0]); niter += 1
            xarr = z1[0]; yarr = z1[1];   # end of sub_size_wenss_getnum
            if noutold == nout: break

          # flag in the 'unresolved' sources. returns flag array, True ==> unresolved
          logsnr=N.log10(snr)
          dumr = N.sqrt(s_c[0]*s_c[0]+s_c[1]*s_c[1]/(snr*snr))
          med = s_dm[0]+s_dm[1]*logsnr+s_dm[2]*(logsnr*logsnr)
          f_sclip[idx] = N.abs((nbeam-med)/(med*dumr)) < N.array([kappa2]*num)
        f_s = f_sclip[0]*f_sclip[1]

        # Add bright sources
        if bright_snr_cut != None:
            if bright_snr_cut < 20.0:
                bright_snr_cut = 20.0
            bright_srcs = N.where(snr >= bright_snr_cut)
            if len(bright_srcs[0]) > 0:
                f_s[bright_srcs] = True

          # now make plots
#           if plot:
#             bb=[b1, b2]
#             pl.subplot(211+idx)
#             pl.semilogx(s1, bb[idx], 'og')
#             f0=f_sclip[idx][index.argsort()]
#             sf=[n for i, n in enumerate(s1) if f0[i]]
#             b1f=[n for i, n in enumerate(bb[idx]) if f0[i]]
#             pl.semilogx(sf, b1f, 'or')
#             pl.semilogx(snr,med,'-')
#             pl.semilogx(snr,med+med*dumr*(N.array([kappa2]*num)),'-')
#             pl.semilogx(snr,med-med*dumr*(N.array([kappa2]*num)),'-')
#             pl.title(' axis ' + str(idx))
#
        return f_s[index.argsort()]

##################################################################################################
    def av_psf(self, g_gauls, beam, flag):
        """ calculate how much the SNR-weighted sizes of unresolved sources differs from the
        synthesized beam. Same as av_psf.f in fBDSM."""
        from math import sqrt

        bmaj = N.asarray(g_gauls[4])
        bmin = N.asarray(g_gauls[5])
        bpa = N.asarray(g_gauls[6])
        wt = N.asarray(g_gauls[1])/N.asarray(g_gauls[8])
        flagwt = wt*flag
        sumwt = N.sum(flagwt)
        w1 = N.sum(flagwt*flagwt)
        wtavbm = N.array([N.sum(bmaj*flagwt), N.sum(bmin*flagwt), N.sum(bpa*flagwt)])/sumwt
        dumrar = N.array([N.sum(bmaj*bmaj*flagwt), N.sum(bmin*bmin*flagwt), N.sum(bpa*bpa*flagwt)])
        dd = sumwt*sumwt-w1
        wtstdbm = N.sqrt((dumrar - wtavbm*wtavbm*sumwt)*sumwt/dd)

        avpa  = N.sum(bpa*flagwt-180.0*flagwt*N.array(bpa >= 90))/sumwt
        stdpa = N.sum(bpa*flagwt+(180.0*180.0-360.0*bpa)*flagwt*N.array(bpa >= 90))
        stdpa = sqrt(abs((stdpa-avpa*avpa*sumwt)*sumwt/dd))
        if stdpa < wtstdbm[2]:
            wtstdbm[2] = stdpa
            wtavbm[2] = avpa

        return (wtavbm - N.array([beam[0]*3600.0, beam[1]*3600.0, beam[2]]))/wtstdbm

##################################################################################################
    def get_voronoi_generators(self, g_gauls, generators, gencode, snrcut, snrtop, snrbot, snrcutstack):
        """This gets the list of all voronoi generators. It is either the centres of the brightest
        sources, or is imported from metadata (in future)."""
        from math import sqrt

        num=len(g_gauls[0])
        snr=N.asarray(g_gauls[1])/N.asarray(g_gauls[8])

        index=snr.argsort()
        snr_incr = snr[index]
        snr = snr_incr[::-1]
        x = N.asarray(g_gauls[2])[index]
        y = N.asarray(g_gauls[3])[index]

        cutoff = 0
        if generators == 'calibrators' or generators == 'field':
            if gencode != 'file':
                gencode = 'list'
            if gencode == 'list':
                cutoff = int(round(num*(snrtop)))
                if cutoff > len(snr):
                    cutoff = len(snr)
                # Make sure we don't fall below snrcutstack (SNR cut for stacking of PSFs), since
                # it makes no sense to make tiles with generators that fall below this cut.
                if snr[cutoff-1] < snrcutstack:
                    cutoff = num - snr_incr.searchsorted(snrcutstack)

        if generators == 'calibrators':
            if gencode == 'file':
                raise NotImplementedError, "gencode=file not yet implemented."

        x1 = x.tolist()
        y1 = y.tolist()
        x1.reverse()
        y1.reverse()
        snr1 = snr.tolist()
        vorogenP = N.asarray([x1[0:cutoff], y1[0:cutoff], snr1[0:cutoff]])

        vorogenS = None

        return vorogenP, vorogenS

##################################################################################################
    def edit_vorogenlist(self, vorogenP, frac):
        """ Edit primary voronoi generator list. Each tile has a tile centre and can
        have more than one generator to be averaged. tile_list is a list of arrays, indexed
        by the tile number and each array is an array of numbers in the ngen list which are
        the generators in that tile. xtile, ytile and snrtile are arrays of length number_of_tiles
        and have x,y,snr of each tile. Group together generators
        if closer than a fraction of dist to third closest."""

        xgen, ygen, snrgen = vorogenP
        flag = N.zeros(len(xgen))
        coord=N.array([xgen,ygen]).transpose()
        tile_list = []
        tile_coord = []; tile_snr = []
        for i in range(len(xgen)):
            dist = N.array(map(lambda t: func.dist_2pt(coord[i], t), coord))
            indi = N.argsort(dist)
            sortdist = dist[indi]
            if sortdist[1] < frac * sortdist[2]:    # first is the element itself
              if flag[indi[1]] + flag[i] == 0:   #  not already deleted from other pair
                tile_list.append([i, indi[1]])
                tile_coord.append((coord[i]*snrgen[i]+coord[indi[1]]*snrgen[indi[1]])/(snrgen[i]+snrgen[indi[1]]))
                tile_snr.append(snrgen[i]+snrgen[indi[1]])
                flag[i] = 1
                flag[indi[1]] = 1
            else:
              if len(dist) > 3:
                if sortdist[1]+sortdist[2] < 2.0*frac*sortdist[3]: # for 3 close-by sources
                  in1=indi[1]
                  in2=indi[2]
                  if flag[in1]+flag[in2]+flag[i] == 0: # not already deleted from others
                    tile_list.append([i, in1, in2])
                    tile_coord.append((coord[i]*snrgen[i]+coord[in1]*snrgen[in1]+coord[in2]*snrgen[in2]) \
                               /(snrgen[i]+snrgen[in1]+snrgen[in2]))
                    tile_snr.append(snrgen[i]+snrgen[in1]+snrgen[in2])
                    flag[i] = 1
                    flag[in1] = 1
                    flag[in2] = 1
              else:
                tile_list.append([i])
                tile_coord.append(coord[i])
                tile_snr.append(snrgen[i])

        # Assign any leftover generators
        for i in range(len(xgen)):
            if flag[i] == 0:
                tile_list.append([i])
                tile_coord.append(coord[i])
                tile_snr.append(snrgen[i])

        return tile_list, tile_coord, tile_snr

##################################################################################################
    def tess_simple(self, vorogenP, wts, tess_sc, tess_fuzzy, shape):
        """ Simple tesselation """

        xgen, ygen, snrgen = vorogenP
        volrank = _pytess.pytess_simple(shape[0], shape[1], xgen, ygen, snrgen, \
                  wts, tess_fuzzy, tess_sc)

        return volrank

##################################################################################################
    def tess_roundness(self, vorogenP, tess_sc, tess_fuzzy, shape):
        """ Tesselation, modified to make the tiles more round. """

        xgen, ygen, snrgen = vorogenP
        volrank = _pytess.pytess_roundness(shape[0], shape[1], xgen, ygen, snrgen, \
                  tess_fuzzy, tess_sc)

        return volrank

##################################################################################################
    def pixintile(self, tilecoord, pixel, tess_method, wts, tess_sc, tess_fuzzy):
        """  This has routines to find out which tile a given pixel belongs to. """

        if tess_method == 'roundness':
          #tilenum = pytess_roundness(tilecoord, pixel, wts, tess_sc, tess_fuzzy)
          print " Not yet implemented !!!! "
          return 0
        else:
          xgen, ygen = tilecoord
          xgen = N.asarray(xgen)
          ygen = N.asarray(ygen)
          ngen = len(xgen)
          i,j = pixel
          dist = N.sqrt((i-xgen)*(i-xgen)+(j-ygen)*(j-ygen))/wts
          minind = dist.argmin()

          if tess_sc == 's':
            tilenum=minind
          else:
            print " Not yet implemented !!!! "

        return tilenum

##################################################################################################
    def tesselate(self, vorogenP, vorogenS, tile_prop, tess_method, tess_sc, tess_fuzzy, generators, gencode, shape):
        """ Various ways of tesselating. If generators='calibrator', no need to tesselate, just get
        modified list based on very nearby sources. If generators='field' then tesselate. The image
        is tesselated based on tile_prop. """

        wtfn={'unity' : lambda x : N.ones(len(x)), \
              'log10' : N.log10, \
              'sqrtlog10' : lambda x : N.sqrt(N.log10(x)), \
              'roundness' : N.array}

        tile_list, tile_coord, tile_snr = tile_prop
        xt = self.trans_gaul(tile_coord)[0]
        yt = self.trans_gaul(tile_coord)[1]
        vorogenT = xt, yt, tile_snr

        wt_fn = wtfn[tess_method]
        wts = wt_fn(tile_snr)

        if tess_method == 'roundness':
            volrank = self.tess_roundness(vorogenT, tess_sc, tess_fuzzy, shape)
        else:
            volrank = self.tess_simple(vorogenT, wts, tess_sc, tess_fuzzy, shape)

        return volrank, wts

##################################################################################################
    def edit_tile(self, ltnum, g_gauls, flag_unresolved, snrcutstack, volrank, tile_prop, tess_sc, \
                  tess_fuzzy, wts, tess_method, plot):
        """ Looks at tiles with no (or one) unresolved source inside it and deletes it and recomputes
           the tiling. For now, does not recompute since we wont use the rank for those pixels anyway."""

        if ltnum > 1: raise NotImplementedError, "NOT YET IMPLEMENTED FOR LTNUM>1"

        tile_list, tile_coord, tile_snr = tile_prop
        tr_gaul = self.trans_gaul(g_gauls)
        tr=[n for i, n in enumerate(tr_gaul) if flag_unresolved[i] and n[1]/n[8] >= snrcutstack]
        ntile = len(tile_list)
        ngenpertile=N.zeros(ntile)
        for itile in range(ntile):
            tile_gauls = [n for n in tr if volrank[int(round(n[2])),int(round(n[3]))]-1 \
                       == itile]
            ngenpertile[itile]=len(tile_gauls)
        new_n = N.sum(ngenpertile >= ltnum)

    # prepare list of good tiles to pass to pixintile
        goodtiles = N.array(N.where(ngenpertile >= ltnum)[0])
        new_n = len(goodtiles)
        tile_coord_n = [n for i,n in enumerate(tile_coord) if i in goodtiles]
        wts_n = [n for i,n in enumerate(wts) if i in goodtiles]

        r2t = N.zeros(ntile, dtype=int)
        entry = -1
        for itile in range(ntile):
          if ngenpertile[itile] >= ltnum:
              r2t[itile] = itile
          else:
              pixel = tile_coord[itile]
              tilenum = self.pixintile(self.trans_gaul(tile_coord_n), pixel, tess_method, wts_n, tess_sc, tess_fuzzy)
              r2t[itile] = tilenum
        for itile in range(new_n):
            num = N.sum(r2t == itile)
            if num == 0:
              minarr = -999
              while minarr != itile:
                arr = N.where(r2t > itile)[0]
                minarr = r2t[arr].min()-1
                for i in arr: r2t[i]=r2t[i]-1

        n_tile_list = []; n_tile_coord = []; n_tile_snr = []
        for itile in range(new_n):
          ind = N.where(r2t == itile)[0]; ind1 = []
          for i in ind: ind1 = ind1 + tile_list[i]
          n_tile_list.append(ind1)
          snrs = N.array([tile_snr[i] for i in ind])
          coords = N.array([tile_coord[i] for i in ind])
          n_tile_snr.append(N.sum(snrs))
          n_tile_coord.append(N.sum([snrs[i]*coords[i] for i in range(len(snrs))], 0)/N.sum(snrs))

        ngenpertile=N.zeros(new_n)
        for itile in range(new_n):
            tile_gauls = [n for n in tr if r2t[volrank[int(round(n[2])),int(round(n[3]))]-1] \
                       == itile]
            ngenpertile[itile]=len(tile_gauls)
        tile_prop = n_tile_list, n_tile_coord, n_tile_snr

        return ngenpertile, tile_prop, r2t

##################################################################################################
    def stackpsf(self, image, beam, g_gauls, wts, cdelt, factor):
        """ Stacks all the images of sources in the gaussian list gauls from image, out to
        a factor times the beam size. Currently the mask is for the whole image but need to
        modify it for masks for each gaussian. These gaussians are supposed to be relatively
        isolated unresolved sources. Cut out an image a big bigger than facXbeam and imageshift
        to nearest half pixel and then add.

        Does not handle masks etc well at all. Masks for image for blanks, masks for \
        islands, etc."""

        gxcens_pix = g_gauls[2]
        gycens_pix = g_gauls[3]
        peak = g_gauls[1]

        psfimsize = int(round(max(beam[0], beam[1])/max(cdelt[0], cdelt[1]) * factor))    # fac X fwhm; fac ~ 2
        psfimage = N.zeros((psfimsize, psfimsize))
        cs2=cutoutsize2 = int(round(psfimsize*(1. + 2./factor)/2.))  # size/2. factor => to avoid edge effects etc
        cc = cutoutcen_ind=[cs2, cs2]
        cpsf=cen_psf_ind = N.array([int(round(psfimsize))/2]*2)
        wt=0.

        num=len(gxcens_pix)
        for isrc in range(num):   #  MASK !!!!!!!!!!!
            wt += wts[isrc]
            gcp=N.array([gxcens_pix[isrc], gycens_pix[isrc]])
            gcen_ind=gcp-1
            rc=rcen_ind = N.asarray(N.round(gcen_ind), dtype=int)
            shift=cc-(gcen_ind-(rc-cs2))
            cutimage = image[rc[0]-cs2:rc[0]+cs2,rc[1]-cs2:rc[1]+cs2]
            if len(cutimage.shape) == 3: cutimage=cutimage[:,:,0]
            if 0 not in cutimage.shape:
                if sum(sum(N.isnan(cutimage))) == 0:
                  im_shift = func.imageshift(cutimage, shift)
                  im_shift = im_shift/peak[isrc]*wts[isrc]
                  subim_shift = im_shift[cc[0]-cpsf[0]:cc[0]-cpsf[0]+psfimsize,cc[1]-cpsf[1]:cc[1]-cpsf[1]+psfimsize]
                  if subim_shift.shape == psfimage.shape:
                    # Check shapes, as they can differ if source is near edge of image.
                    # If they do differ, don't use that source (may be distorted).
                    psfimage += subim_shift
        psfimage = psfimage/wt

        return psfimage

##################################################################################################
    def psf_in_tile(self, image, beam, g_gauls, cdelt, factor, snrcutstack, volrank, \
                    tile_prop, plot, img):
        """ For each tile given by tile_prop, make a list of all gaussians in the constituent tesselations
        and pass it to stackpsf with a weight for each gaussian, to calculate the average psf per tile.

        Should define weights inside a tile to include closure errors """
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Psf_Vary")

        tile_list, tile_coord, tile_snr = tile_prop
        tr_gaul = self.trans_gaul(g_gauls)
        tr=[n for i, n in enumerate(tr_gaul)]# if n[1]/n[8] >= snrcutstack]
        ntile = len(tile_list)
        psfimages = []
        psfcoords = []
        psfratio = [] # ratio of peak flux to total flux
        psfratio_aper = [] # ratio of peak flux to aperture flux
        srcpertile = N.zeros(ntile)
        snrpertile = N.zeros(ntile)
        xt, yt = N.transpose(tile_coord)

        if plot:
          pl.figure(None)
          colours=['b','g','r','c','m','y','k']*(len(xt)/7+1)
          pl.axis([0.0, image.shape[0], 0.0, image.shape[1]])
          pl.title('Tesselated image with tile centres and unresolved sources')
          for i in range(ntile):
            pl.plot([xt[i]], [yt[i]], 'D'+colours[i])
            pl.text(xt[i], yt[i], str(i))

        for itile in range(ntile):
            tile_gauls = [n for n in tr if volrank[int(round(n[2])),int(round(n[3]))]-1 \
                       == itile]
            t_gauls = self.trans_gaul(tile_gauls)
            srcpertile[itile] = len(tile_gauls)
            if plot:
              pl.plot(t_gauls[2], t_gauls[3], 'x'+'k', mew=1.3)#colours[itile])
              for i, ig in enumerate(t_gauls[2]):
                xx=[xt[itile], ig]
                yy=[yt[itile], t_gauls[3][i]]
                pl.plot(xx,yy,'-'+colours[itile])
            wts = N.asarray(t_gauls[1])/N.asarray(t_gauls[8])             # wt is SNR
            snrpertile[itile] = sum(wts)
            mylog.info('PSF tile #%i (center = %i, %i): %i unresolved sources, SNR = %.1f' %
                (itile, xt[itile], yt[itile], srcpertile[itile], snrpertile[itile]))
            a = self.stackpsf(image, beam, t_gauls, wts, cdelt, factor)
            psfimages.append(a)
            psfcoords.append([sum(N.asarray(t_gauls[2])*wts)/sum(wts), sum(N.asarray(t_gauls[3])*wts)/sum(wts)])

            # Find peak/total flux ratio for sources in tile. If an aperture is given,
            # use the aperture flux as well.
            # t_gauls[0] is source_id
            src_ratio = []
            src_wts = []
            src_ratio_aper = []
            src_wts_aper = []
            for gt in tile_gauls:
                src = img.sources[gt[0]]
                if img.aperture != None:
                    src_ratio_aper.append(src.peak_flux_max / src.aperture_flux)
                    src_wts_aper.append(src.total_flux / src.aperture_fluxE)
                src_ratio.append(src.peak_flux_max / src.total_flux)
                src_wts.append(src.total_flux / src.total_fluxE)
            if img.aperture != None:
                psfratio_aper.append(sum(N.asarray(src_ratio_aper)*src_wts_aper)/sum(src_wts_aper))
            else:
                psfratio_aper.append(0.0)
            psfratio.append(sum(N.asarray(src_ratio)*src_wts)/sum(src_wts))

        totpsfimage = psfimages[0]*snrpertile[0]
        for itile in range(1,ntile):
            totpsfimage += psfimages[itile]*snrpertile[itile]
        totpsfimage = totpsfimage/sum(snrpertile)

        if plot:
          pl.imshow(N.transpose(volrank), origin='lower', interpolation='nearest'); pl.colorbar()

        if plot:
         pl.figure(None)
         pl.clf()
         ax = pl.subplot(1,1,1)
         pax = ax.get_position()
         start = N.array((pax.xmin, pax.ymin))
         stop = N.array((pax.xmax, pax.ymax))
         plaxis = pl.axis([0, image.shape[0], 0, image.shape[1]])
         pl.title('Stacked psf for each tile')
         for itile in range(ntile):
           im=psfimages[itile]
           sz=0.07
           spt = int(round(snrpertile[itile]*10))/10.
           titl='n='+str(int(round(srcpertile[itile])))+'; SNR='+str(spt)
           posn=[psfcoords[itile][0], psfcoords[itile][1]]
           normposn=N.array(stop-start, dtype=float)/N.array(image.shape[0:2])*posn+start
           a=pl.axes([normposn[0]-sz/2., normposn[1]-sz/2., sz, sz])
           pl.contour(im,15)
           pl.title(titl, fontsize='small')
           pl.setp(a, xticks=[], yticks=[])
         pl.show()

        return psfimages, psfcoords, totpsfimage, psfratio, psfratio_aper


##################################################################################################
    def interp_shapcoefs(self, nmax, tr_psf_cf, psfcoords, imshape, compress, plot):
        """Interpolate using natgrid.

        Check to see if variation is significant.
        """
        x, y = N.transpose(psfcoords)
        index = [(i,j) for i in range(nmax+1) for j in range(nmax+1-i)]
        xi=x
        yi=y
        xo=N.arange(0.0,round(imshape[0]), round(compress))
        yo=N.arange(0.0,round(imshape[1]), round(compress))
        rgrid=nat.Natgrid(xi,yi,xo,yo)
        p={}
        for coord in index:
            z = N.array(tr_psf_cf[coord])    # else natgrid cant deal with noncontiguous memory
            p[coord] = rgrid.rgrd(z)

#         if plot:
#           for i,coord in enumerate(index):
#             if i % 36 == 0:
#               pl.figure(None)
#               pl.clf()
#               title = 'Interpolated shapelet coefficients'
#               if i>0: title = title+' (cont)'
#               pl.suptitle(title)
#             pl.subplot(6,6,(i%36)+1)
#             pl.title(str(coord))
#             pl.plot(xi/compress, yi/compress, 'xk')
#             pl.imshow(p[coord], interpolation='nearest')
#             pl.colorbar()

        return p, xo, yo

##################################################################################################
    def interp_prop(self, prop, psfcoords, imshape, compress=1):
        """Interpolate using natgrid.

        Should check to see if variation is significant.
        """
        x, y = N.transpose(psfcoords)
        xi=x
        yi=y
        xo=N.arange(0.0,round(imshape[0]), round(compress))
        yo=N.arange(0.0,round(imshape[1]), round(compress))
        rgrid=nat.Natgrid(xi,yi,xo,yo)
        prop_int = rgrid.rgrd(prop)
        return prop_int

##################################################################################################
    def create_psf_grid(self, psf_coeff_interp, imshape, xgrid, ygrid, skip, nmax, psfshape, basis, beta,
        cen, totpsfimage, plot):
        """ Creates a image with the gridded interpolated psfs. xgrid and ygrid are 1d numpy arrays
        with the x and y coordinates of the grids. """

#         if plot:
#           plnum=N.zeros(2)
#           for i in range(2):
#             dum=pl.figure(None)
#             plnum[i]=dum.number
#             pl.clf()
#             if i == 0: pl.suptitle('Gridded psfs')
#             if i == 1: pl.suptitle('Gridded residual psfs')
#             ax = pl.subplot(1,1,1)
#             plaxis = pl.axis([0, imshape[0], 0, imshape[1]])
#             pax = ax.get_position()
#             start = N.array((pax.xmin, pax.ymin))
#             stop = N.array((pax.xmax, pax.ymax))
#             sz=0.07
        mask=N.zeros(psfshape, dtype=bool)   # right now doesnt matter
        xg=xgrid[::skip+1]
        yg=ygrid[::skip+1]
        index = [(i,j) for i in range(0,len(xgrid),skip+1) for j in range(0,len(ygrid),skip+1)]
        xy = [(i,j) for i in xgrid[::skip+1] for j in ygrid[::skip+1]]
        blah=[]
        for i, coord in enumerate(index):
            maxpsfshape = [0, 0]
            for k in psf_coeff_interp:
                if k[0]+1 > maxpsfshape[0]:
                    maxpsfshape[0] = k[0]+1
                if k[1]+1 > maxpsfshape[1]:
                    maxpsfshape[1] = k[1]+1
            cf = N.zeros(maxpsfshape)
            for k in psf_coeff_interp:
                cf[k]=psf_coeff_interp[k][coord]
            cf = N.transpose(cf)
            psfgridim = sh.reconstruct_shapelets(psfshape, mask, basis, beta, cen, nmax, cf)
            blah.append(psfgridim)

#             if plot:
#               for j in range(2):
#                 pl.figure(plnum[j])
#                 posn = [xy[i][0], xy[i][1]]
#                 normposn =N.array(stop-start, dtype=float)/N.array(imshape[0:2])*posn+start
#                 a=pl.axes([normposn[0]-sz/2., normposn[1]-sz/2., sz, sz])
#                 if j == 0: pl.contour(psfgridim,15)
#                 if j == 1: pl.contour(psfgridim-totpsfimage,15)
#                 pl.setp(a, xticks=[], yticks=[])
#                 pl.colorbar()
#         if plot:
#           pl.figure(plnum[0])
#           pl.figure(plnum[1])
#
        return blah

##################################################################################################
    def blur_image(self, im, n, ny=None) :
        """ blurs the image by convolving with a gaussian kernel of typical
            size n. The optional keyword argument ny allows for a different
            size in the y direction.
        """
        from scipy.ndimage import gaussian_filter

        sx = n
        if ny != None:
            sy = ny
        else:
            sy = n
        improc = gaussian_filter(im, [sy, sx])
        return improc
