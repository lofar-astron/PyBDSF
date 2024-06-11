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
from __future__ import print_function
from __future__ import absolute_import

from .image import *
from . import mylogger
from . import statusbar
from . import has_pl
if has_pl:
    import matplotlib.pyplot as pl
import scipy.ndimage as nd
from . import multi_proc as mp
import itertools


class Op_gausfit(Op):
    """Fit a number of 2D gaussians to each island.

    The results of the fitting are stored in the Island
    structure itself as a list of Gaussian objects (gaul) and a
    list of flagged gaussians (fgaul).

    Prerequisites: module islands should be run first.
    """
    def __call__(self, img):
        from . import functions as func

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
        if not opts.quiet and not opts.verbose_fitting:
            bar.start()
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
        img_simple.beam2pix = img.beam2pix
        img_simple.beam = img.beam

        # Next, define the weights to use when distributing islands among cores.
        # The weight should scale with the processing time. At the moment
        # we use the island area, but other parameters may be better.
        weights = []
        for isl in img.islands:
            weights.append(isl.size_active)

        # Now call the parallel mapping function. Returns a list of
        # [gaul, fgaul] for each island.  If ncores is 1, use the
        # standard Python map function -- this helps with debugging in
        # some circumstances
        if opts.ncores == 1:
            gaus_list = map(func.eval_func_tuple,
                            zip(itertools.repeat(self.process_island),
                                img.islands, itertools.repeat(img_simple),
                                itertools.repeat(opts)))
        else:
            gaus_list = mp.parallel_map(func.eval_func_tuple,
                                        zip(itertools.repeat(self.process_island),
                                            img.islands, itertools.repeat(img_simple),
                                            itertools.repeat(opts)),
                                        numcores=opts.ncores, bar=bar, weights=weights)
        gaus_list = list(gaus_list)

        for isl in img.islands:
            # Now convert gaussians into Gaussian objects and store
            idx = isl.island_id
            gaul = gaus_list[idx][0]
            fgaul = gaus_list[idx][1]
            dgaul = []
            if len(gaul) > 0:
                gidx = gaul[-1][0]  # save last index value for use with fgaul below
            else:
                gidx = 0
            gaul = [Gaussian(img, par, idx, gidx)
                    for (gidx, par) in enumerate(gaul)]

            if len(gaul) == 0:
                # No good Gaussians were fit. In this case, make a dummy
                # Gaussian located at the island center so
                # that the source may still be included in output catalogs.
                # These dummy Gaussians all have an ID of -1. They do not
                # appear in any of the source or island Gaussian lists except
                # the island dgaul list.
                if opts.src_ra_dec is not None:
                    # Center the dummy Gaussian on the user-specified source position
                    posn_isl = (int(isl.shape[0]/2.0), int(isl.shape[1]/2.0))
                    posn_img = (int(isl.shape[0]/2.0 + isl.origin[0]), int(isl.shape[1]/2.0 + isl.origin[1]))
                    par = [isl.image[posn_isl], posn_img[0], posn_img[1], 0.0, 0.0, 0.0]
                else:
                    # Center the dummy Gaussian on the maximum pixel
                    posn = N.unravel_index(N.argmax(isl.image*~isl.mask_active), isl.shape) + N.array(isl.origin)
                    par = [isl.max_value, posn[0], posn[1], 0.0, 0.0, 0.0]
                dgaul = [Gaussian(img, par, idx, -1)]

            # Now make the list of flagged Gaussians, if any
            fgaul = [Gaussian(img, par, idx, gidx + gidx2 + 1, flag)
                     for (gidx2, (flag, par)) in enumerate(fgaul)]

            isl.gaul = gaul
            isl.fgaul = fgaul
            isl.dgaul = dgaul

        gaussian_list = [g for isl in img.islands for g in isl.gaul]
        img.gaussians = gaussian_list

        # Put in the serial number of the gaussians for the whole image
        n = 0
        nn = 0
        tot_flux = 0.0
        for isl in img.islands:
            m = 0
            for g in isl.gaul:
                n += 1
                m += 1
                g.gaus_num = n - 1
                tot_flux += g.total_flux
            for dg in isl.dgaul:
                nn -= 1
                dg.gaus_num = nn

            isl.ngaus = m
        img.ngaus = n
        img.total_flux_gaus = tot_flux

        mylogger.userinfo(mylog, "Total number of Gaussians fit to image",
                          str(n))
        if not img._pi and not img.waveletimage:
            mylogger.userinfo(mylog, "Total flux density in model", '%.3f Jy' %
                              tot_flux)

        # Check if model flux is very different from sum of flux in image
        if img.ch0_sum_jy > 0 and not img._pi:
            if img.total_flux_gaus/img.ch0_sum_jy < 0.5 or \
                    img.total_flux_gaus/img.ch0_sum_jy > 2.0:
                mylog.warn('Total flux density in model is %0.2f times sum of pixels '
                           'in input image. Large residuals may remain.' %
                           (img.total_flux_gaus/img.ch0_sum_jy,))

        # Check if there are many Gaussians with deconvolved size of 0 in one
        # axis but not in the other. Don't bother to do this for wavelet images.
        fraction_1d = self.check_for_1d_gaussians(img)
        if fraction_1d > 0.5 and img.beam is not None and not img.waveletimage:
            mylog.warn("After deconvolution, more than 50% of Gaussians are "
                       "1-D. Unless you're fitting an extended source, "
                       "beam may be incorrect.")

        img.completed_Ops.append('gausfit')
        return img

    def process_island(self, isl, img, opts=None):
        """Processes a single island.

        Returns a list of the best-fit Gaussians and flagged Gaussians.
        """
        from . import functions as func

        if opts is None:
            opts = img.opts
        iter_ngmax = 10
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

        size = isl.size_active/img.pixel_beamarea()*2.0   # 2.0 roughly corrects for thresh_isl
        if opts.verbose_fitting:
            print("Fitting isl #", isl.island_id, '; # pix = ', N.sum(~isl.mask_active), '; size = ', size)

        if size > maxsize:
            tosplit = func.isl_tosplit(isl, opts)
            if opts.split_isl and tosplit[0] > 0:
                n_subisl, sub_labels = tosplit[1], tosplit[2]
                gaul = []
                fgaul = []
                if opts.verbose_fitting:
                    print('SPLITTING ISLAND INTO ', n_subisl, ' PARTS FOR ISLAND ', isl.island_id)
                for i_sub in range(n_subisl):
                    islcp = isl.copy(img.pixel_beamarea())
                    islcp.mask_active = N.where(sub_labels == i_sub+1, False, True)
                    islcp.mask_noisy = N.where(sub_labels == i_sub+1, False, True)
                    size_subisl = (~islcp.mask_active).sum()/img.pixel_beamarea()*2.0
                    if opts.peak_fit and size_subisl > peak_size:
                        sgaul, sfgaul = self.fit_island_iteratively(img, islcp, iter_ngmax=iter_ngmax, opts=opts)
                    else:
                        sgaul, sfgaul = self.fit_island(islcp, opts, img)
                    gaul = gaul + sgaul
                    fgaul = fgaul + sfgaul
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

        Note: "fitok" indicates whether fit converged
               and one or more flagged Gaussians indicate
               that significant residuals remain (peak > thr).
        """
        from ._cbdsm import MGFunction
        from . import functions as func
        from .const import fwsig

        verbose = opts.verbose_fitting
        if verbose:
            print('Entering fit_island in verbose mode')

        if ffimg is None:
            fit_image = isl.image-isl.islmean
        else:
            fit_image = isl.image-isl.islmean-ffimg

        fcn = MGFunction(fit_image, isl.mask_active, 1)
        # For fitting, use img.beam instead of img.pixel_beam, as we want
        # to pick up the wavelet beam (img.pixel_beam is not changed for
        # wavelet images, but img.beam is)
        beam = N.array(img.beam2pix(img.beam))
        beam = (beam[0]/fwsig, beam[1]/fwsig, beam[2]+90.0)  # change angle from +y-axis to +x-axis and FWHM to sigma

        if abs(beam[0]/beam[1]) < 1.1:
            beam = (1.1*beam[0], beam[1], beam[2])

        thr1 = isl.mean + opts.thresh_isl*isl.rms
        thr0 = thr1
        g3_only = opts.fix_to_beam
        peak = fcn.find_peak()[0]
        dof = isl.size_active
        shape = isl.shape
        size = isl.size_active/img.pixel_beamarea()*2.0
        gaul = []
        iter = 0
        ng1 = 0
        if ini_gausfit is None:
            ini_gausfit = opts.ini_gausfit

        if ini_gausfit not in ['default', 'simple', 'nobeam']:
            ini_gausfit = 'default'
        if ini_gausfit == 'simple' and ngmax is None:
            ngmax = 25
        if ini_gausfit == 'default' or opts.fix_to_beam:
            gaul, ng1, ngmax = self.inigaus_fbdsm(isl, thr0, beam, img)
            if len(gaul) > 25:
                ini_gausfit = 'simple'
                gaul = []
                ng1 = 0
                ngmax = 25
        if ini_gausfit == 'nobeam' and not opts.fix_to_beam:
            gaul = self.inigaus_nobeam(isl, thr0, beam, img)
            ng1 = len(gaul)
            ngmax = ng1+2
        if verbose:
            print('Initializing, ini_gausfit is', ini_gausfit, 'gaul =', gaul, 'ngmax =', ngmax)
        while iter < 5:
            iter += 1
            if verbose:
                print('In Gaussian flag loop, iter =', iter)
            fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, ini_gausfit, ngmax, verbose, g3_only)
            if verbose:
                print('Calling flag_gaussians')
            gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                              beam, thr0, peak, shape, isl.mask_active,
                                              isl.image, size)
            if verbose:
                print('Leaving flag_gaussians')
            ng1 = len(gaul)
            if fitok and len(fgaul) == 0:
                break
        if (not fitok or len(gaul) == 0) and ini_gausfit != 'simple':
            if verbose:
                print('Using simple method instead')
            # If fits using default or nobeam methods did not work,
            # try using simple instead
            gaul = []
            iter = 0
            ng1 = 0
            ngmax = 25
            while iter < 5:
                iter += 1
                fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, 'simple', ngmax, verbose, g3_only)
                gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                                  beam, thr0, peak, shape, isl.mask_active,
                                                  isl.image, size)
                ng1 = len(gaul)
                if fitok and len(fgaul) == 0:
                    break
        sm_isl = nd.binary_dilation(isl.mask_active)
        if (not fitok or len(gaul) == 0) and N.sum(~sm_isl) >= img.minpix_isl:
            if verbose:
                print('Fit still not OK, shrinking')
            # If fitting still fails, shrink the island a little and try again
            fcn = MGFunction(fit_image, nd.binary_dilation(isl.mask_active), 1)
            gaul = []
            iter = 0
            ng1 = 0
            ngmax = 25
            while iter < 5:
                iter += 1
                fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, 'simple', ngmax, verbose, g3_only)
                gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                                  beam, thr0, peak, shape, isl.mask_active,
                                                  isl.image, size)
                ng1 = len(gaul)
                if fitok and len(fgaul) == 0:
                    break
        lg_isl = nd.binary_erosion(isl.mask_active)
        if (not fitok or len(gaul) == 0) and N.sum(~lg_isl) >= img.minpix_isl:
            if verbose:
                print('Fit still not OK, expanding')
            # If fitting still fails, expand the island a little and try again
            fcn = MGFunction(fit_image, nd.binary_erosion(isl.mask_active), 1)
            gaul = []
            iter = 0
            ng1 = 0
            ngmax = 25
            while iter < 5:
                iter += 1
                fitok = self.fit_iter(gaul, ng1, fcn, dof, beam, thr0, iter, 'simple', ngmax, verbose, g3_only)
                gaul, fgaul = self.flag_gaussians(fcn.parameters, opts,
                                                  beam, thr0, peak, shape, isl.mask_active,
                                                  isl.image, size)
                ng1 = len(gaul)
                if fitok and len(fgaul) == 0:
                    break

        if not fitok or len(gaul) == 0:
            # If all else fails, try to use moment analysis
            if verbose:
                print('All else has failed, trying moment analysis')
            inisl = N.where(~isl.mask_active)
            mask_id = N.zeros(isl.image.shape, dtype=N.int32) - 1
            mask_id[inisl] = isl.island_id
            try:
                pixel_beamarea = img.pixel_beamarea()
                mompara = func.momanalmask_gaus(fit_image, mask_id, isl.island_id, pixel_beamarea, True)
                mompara[5] += 90.0
                if not N.isnan(mompara[1]) and not N.isnan(mompara[2]):
                    x1 = int(N.floor(mompara[1]))
                    y1 = int(N.floor(mompara[2]))
                    t = (mompara[1]-x1)/(x1+1-x1)
                    u = (mompara[2]-y1)/(y1+1-y1)
                    s_peak = ((1.0-t) * (1.0-u) * fit_image[x1, y1] + t * (1.0-u) * fit_image[x1+1, y1] +
                              t * u * fit_image[x1+1, y1+1] + (1.0-t) * u * fit_image[x1, y1+1])
                    mompara[0] = s_peak
                    par = [mompara.tolist()]
                    par[3] /= fwsig
                    par[4] /= fwsig
                    gaul, fgaul = self.flag_gaussians(par, opts,
                                                      beam, thr0, peak, shape, isl.mask_active,
                                                      isl.image, size)
            except:
                pass

        # Return whatever we got
        if verbose:
            print('Preparing to return')
        isl.mg_fcn = fcn
        gaul = [self.fixup_gaussian(isl, g) for g in gaul]
        fgaul = [(flag, self.fixup_gaussian(isl, g)) for flag, g in fgaul]

        if verbose:
            print('Number of good Gaussians: %i' % (len(gaul),))
            print('Number of flagged Gaussians: %i' % (len(fgaul),))
        return gaul, fgaul

    def fit_island_iteratively(self, img, isl, iter_ngmax=5, opts=None):
        """Fits an island iteratively.

        For large islands, which can require many Gaussians to fit well,
        it is much faster to fit a small number of Gaussians simultaneously
        and iterate. However, this does usually result in larger residuals.
        """
        from . import functions as func
        sgaul = []
        sfgaul = []
        gaul = []
        fgaul = []
        if opts is None:
            opts = img.opts
        thresh_isl = opts.thresh_isl
        thresh = opts.fittedimage_clip
        thr = isl.mean + thresh_isl * isl.rms

        if opts.verbose_fitting:
            print('Iteratively fitting island ', isl.island_id)
        gaul = []
        fgaul = []
        ffimg_tot = N.zeros(isl.shape, dtype=N.float32)
        peak_val = N.max(isl.image - isl.islmean)
        count = 0
        while peak_val >= thr:
            count += 1
            if opts.verbose_fitting:
                print('Iteration %i' % count)
            sgaul, sfgaul = self.fit_island(isl, opts, img, ffimg=ffimg_tot, ngmax=iter_ngmax, ini_gausfit='simple')
            gaul = gaul + sgaul
            fgaul = fgaul + sfgaul

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
                peak_val_prev = peak_val
                peak_val = N.max(isl.image - isl.islmean - ffimg_tot)
                if func.approx_equal(peak_val, peak_val_prev):
                    break
            else:
                break

        if len(gaul) == 0:
            if opts.verbose_fitting:
                # Fitting iteratively did not work -- try normal fit
                print('Iterative fitting failed for', isl.island_id)
            gaul, fgaul = self.fit_island(isl, opts, img, ini_gausfit='default')
        else:
            if opts.verbose_fitting:
                print('Iterative fitting succeeded for', isl.island_id)

        return gaul, fgaul

    def inigaus_fbdsm(self, isl, thr, beam, img):
        """ initial guess for gaussians like in fbdsm """
        from math import sqrt
        from .const import fwsig
        from . import functions as func

        im = isl.image-isl.islmean
        if img.opts.ini_method == 'curvature':
            im_pos = -1.0 * func.make_curvature_map(isl.image-isl.islmean)
            thr_pos = 0.0
        else:
            im_pos = im
            thr_pos = thr
        mask = isl.mask_active
        av = img.clipped_mean
        inipeak, iniposn, im1 = func.get_maxima(im, mask, thr_pos, isl.shape, beam, im_pos=im_pos)
        if len(inipeak) == 0:
            av, stdnew, maxv, maxp, minv, minp = func.arrstatmask(im, mask)
            inipeak = [maxv]
            iniposn = [maxp]
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
                    iniposn.append(N.array(maxp))
                    inipeak.append(maxv)
                    im1 = func.mclean(im1, maxp, beam)

        inipeak = N.array(inipeak)
        iniposn = N.array(iniposn)
        ind = list(N.argsort(inipeak))
        ind.reverse()
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
        from .const import fwsig
        import scipy.ndimage as nd
        from . import functions as func

        im = isl.image-isl.islmean
        if img.opts.ini_method == 'curvature':
            im_pos = -1.0 * func.make_curvature_map(isl.image-isl.islmean)
            thr_pos = 0.0
        else:
            im_pos = im
            thr_pos = -1e9
        mask = isl.mask_active
        av = img.clipped_mean
        inipeak, iniposn, im1 = func.get_maxima(im, mask, thr_pos, isl.shape, beam, im_pos=im_pos)
        npeak = len(iniposn)
        gaul = []

        av, stdnew, maxv, maxp, minv, minp = func.arrstatmask(im, mask)
        mom = func.momanalmask_gaus(isl.image-isl.islmean, isl.mask_active, 0, 1.0, True)
        if npeak <= 1:
            g = (float(maxv), int(round(mom[1])), int(round(mom[2])), mom[3]/fwsig,
                 mom[4]/fwsig, mom[5])
            gaul.append(g)

        if npeak > 1:  # markers start from 1=background, watershed starts from 1=background
            watershed, markers = func.watershed(im, mask=isl.mask_active)
            nshed = N.max(markers)-1  # excluding background
            xm, ym = N.transpose([N.where(markers == i) for i in range(1, nshed+2)])[0]
            coords = [c for c in N.transpose([xm, ym])[1:]]
            alldists = [func.dist_2pt(c1, c2) for c1 in coords for c2 in coords if N.any(c1 != c2)]  # has double
            meandist = N.mean(alldists)  # mean dist between all pairs of markers
            # Find at least some 'compact' sources
            cscale = 3.0
            while True:
                compact = []
                invmask = []
                for ished in range(nshed):
                    shedmask = N.where(watershed == ished+2, False, True) + isl.mask_active  # good unmasked pixels = 0
                    imm = nd.binary_dilation(~shedmask, N.ones((3, 3), int))
                    xbad, ybad = N.where((imm == 1)*(im > im[xm[ished+1], ym[ished+1]]))
                    imm[xbad, ybad] = 0
                    invmask.append(imm)
                    x, y = N.where(imm)
                    xcen, ycen = N.mean(x), N.mean(y)  # good pixels are now = 1
                    dist = func.dist_2pt([xcen, ycen], [xm[ished+1], ym[ished+1]])
                    if dist < max(cscale, meandist/4.0):
                        compact.append(True)  # if not compact, break source + diffuse
                    else:
                        compact.append(False)
                if N.any(compact):
                    break
                else:
                    # Rescale to search for more compact sources
                    cscale *= 1.5

            if not N.all(compact):
                o_avsize = []
                ind = N.where(compact)[0]
                for i in ind:
                    o_avsize.append(N.sum(invmask[i]))
                avsize = sqrt(N.mean(N.array(o_avsize)))
                for i in range(len(compact)):
                    if not compact[i]:  # make them all compact
                        newmask = N.zeros(imm.shape, bool)
                        newmask[max(0, int(xm[i+1]-avsize/2)):min(im.shape[0], int(xm[i+1]+avsize/2)),
                                max(0, int(ym[i+1]-avsize/2)):min(im.shape[1], int(ym[i+1]+avsize/2))] = True
                        invmask[i] = invmask[i]*newmask
            resid = N.zeros(im.shape, dtype=N.float32)  # approx fit all compact ones
            for i in range(nshed):
                size = sqrt(N.sum(invmask))/fwsig
                xf, yf = coords[i][0], coords[i][1]
                p_ini = [im[xf, yf], xf, yf, size, size, 0.0]
                x, y = N.indices(im.shape)
                p, success = func.fit_gaus2d(im*invmask[i], p_ini, x, y)
                resid = resid + func.gaus_2d(p, x, y)
                gaul.append(p)
            resid = im - resid
            if not N.all(compact):  # just add one gaussian to fit whole unmasked island
                maxv = N.max(resid)  # assuming resid has only diffuse emission. can be false
                x, y = N.where(~isl.mask_active)
                xcen = N.mean(x)
                ycen = N.mean(y)
                invm = ~isl.mask_active
                mom = func.momanalmask_gaus(invm, N.zeros(invm.shape, dtype=N.int16), 0, 1.0, True)
                g = (maxv, xcen, ycen, mom[3]/fwsig, mom[4]/fwsig, mom[5]-90.)
                gaul.append(g)
                coords.append([xcen, ycen])

        return gaul

    def fit_iter(self, gaul, ng1, fcn, dof, beam, thr, iter, inifit, ngmax, verbose=1, g3_only=False):
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
        from ._cbdsm import lmder_fit

        if verbose:
            print('Greetings from fit_iter')
        fit = lmder_fit
        beam = list(beam)

        # First drop-in initial gaussians
        # No error-checking here, they MUST fit
        fcn.reset()
        for ig in range(ng1):
            g = gaul[ig]
            self.add_gaussian(fcn, g, dof, g3_only)

        # Do a round of fitting if any initials were provided
        if verbose:
            print('About to call C++ wrapper')
        fitok = True
        if len(gaul) != 0:
            fitok = fit(fcn, final=0, verbose=verbose)

        if verbose:
            print('Returned from the fit')
        # Iteratively add gaussians while there are high peaks
        # in the image and fitting converges
        while fitok:
            peak, coords = fcn.find_peak()
            if peak < thr:  # no good peaks left
                break
            if len(fcn.parameters) < ngmax and iter == 1 and inifit == 'default' and len(gaul) >= ng1+1:
                ng1 = ng1 + 1
                g = gaul[ng1-1]
            else:
                if len(fcn.parameters) < ngmax:
                    g = [peak, coords[0], coords[1]] + beam
                else:
                    break
            fitok &= self.add_gaussian(fcn, g, dof, g3_only)

            fitok &= fit(fcn, final=0, verbose=verbose)

        # And one last fit with higher precision
        # make sure we return False when fitok==False due to lack
        # of free parameters
        fitok &= fit(fcn, final=1, verbose=verbose)

        return fitok

    def add_gaussian(self, fcn, parameters, dof, g3_only=False):
        """Try adding one more gaussian to fcn object.
        It's trying to reduce number of fitted parameters if
        there is not enough DoF left.

        Note: g1 fits amplitude only
              g3 fits amplitude and position
              g6 fits all parameters

        Parameters:
        fcn: MGFunction object
        parameters: initial values for gaussian parameters
        dof: total possible number of fitted parameters
        """
        from ._cbdsm import Gtype

        if g3_only:
            gtype = (Gtype.g3 if fcn.fitted_parameters() + 3 <= dof else None)
        else:
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
        bad = []
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
        from math import sqrt, sin, cos, pi
        from . import functions as func
        import scipy.ndimage as nd

        A, x1, x2, s1, s2, th = g
        s1, s2 = map(abs, [s1, s2])
        flag = 0
        if N.any(N.isnan(g)) or s1 == 0.0 or s2 == 0.0:
            return -1

        if s1 < s2:   # s1 etc are sigma
            ss1 = s2
            ss2 = s1
            th1 = divmod(th+90.0, 180)[1]
        else:
            ss1 = s1
            ss2 = s2
            th1 = divmod(th, 180)[1]
        th1 = th1/180.0*pi
        if ss1 > 1e4 and ss2 > 1e4:
            xbox = 1e9
            ybox = 1e9
        else:
            xbox = 2.0 * (abs(ss1 * cos(th1) * cos(th1)) + abs(ss2 * ss2 / ss1 * sin(th1) * sin(th1))) / \
                   sqrt(cos(th1) * cos(th1) + ss2 * ss2 / ss1 / ss1 * sin(th1) * sin(th1))
            ybox = 2.0 * (abs(ss1 * sin(th1) * sin(th1)) + abs(ss2 * ss2 / ss1 * cos(th1) * cos(th1))) / \
                   sqrt(sin(th1) * sin(th1) + ss2 * ss2 / ss1 / ss1 * cos(th1) * cos(th1))

        # Now check all conditions
        border = opts.flag_bordersize
        x1ok = True
        x2ok = True
        flagmax = False
        if A < opts.flag_minsnr*thr:
            flag += 1
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
            if borx1_1 < 0:
                borx1_1 = 0
            borx1_2 = x1 + border + 1
            if borx1_2 > shape[0]:
                borx1_2 = shape[0]
            if N.any(mask[int(borx1_1):int(borx1_2), int(x2)]):
                flag += 4
            borx2_1 = x2 - border
            if borx2_1 < 0:
                borx2_1 = 0
            borx2_2 = x2 + border + 1
            if borx2_2 > shape[1]:
                borx2_2 = shape[1]
            if N.any(mask[int(x1), int(borx2_1):int(borx2_2)]):
                flag += 8
        if xbox > opts.flag_maxsize_isl*shape[0]:
            flag += 16
        if ybox > opts.flag_maxsize_isl*shape[1]:
            flag += 32
        if s1*s2 > opts.flag_maxsize_bm*beam[0]*beam[1]:
            flag += 64
        if opts.flag_smallsrc:
            if s1*s2 < opts.flag_minsize_bm*beam[0]*beam[1]:
                flag += 128
        if not opts.flag_smallsrc:
            if s1*s2 == 0.:
                flag += 128

        if ss1/ss2 > 2.0:
            # Only check for fairly elliptical Gaussians, as this condition
            # is unreliable for more circular ones.
            ellx, elly = func.drawellipse([A, x1, x2, s1*opts.flag_maxsize_fwhm,
                                           s2*opts.flag_maxsize_fwhm, th])
            pt1 = [N.min(ellx), elly[N.argmin(ellx)]]
            pt2 = [ellx[N.argmax(elly)], N.max(elly)]
            pt3 = [N.max(ellx), elly[N.argmax(ellx)]]
            pt4 = [ellx[N.argmin(elly)], N.min(elly)]
            extremes = [pt1, pt2, pt3, pt4]
            for pt in extremes:
                if N.any(N.isnan(pt)):
                    flag += 256
                    break
                elif pt[0] < 0 or pt[0] >= shape[0] or pt[1] < 0 or pt[1] >= shape[1]:
                    flag += 256
                    break
                elif mask[int(pt[0]), int(pt[1])]:
                    flag += 256
                    break
        return flag

    def fixup_gaussian(self, isl, gaussian):
        """Normalize parameters by adjusting them to the
        proper image coordinates and ensuring that all of
        the implicit conventions (such as bmaj >= bmin) are met.
        """
        np = list(gaussian)

        # Update to the image coordinates
        np[1] += isl.origin[0]
        np[2] += isl.origin[1]

        # Shape values should be positive
        np[3] = abs(np[3])
        np[4] = abs(np[4])

        # First extent is major
        if np[3] < np[4]:
            np[3:5] = np[4:2:-1]
            np[5] += 90

        # Clip position angle
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


class Gaussian(object):
    """Instances of this class are used to store information about
    extracted gaussians in a structured way.
    """

    def __init__(self, img, gaussian, isl_idx, g_idx, flg=0):
        """Initialize Gaussian object from fitting data

        Parameters:
        img: PyBDSM image object
        gaussian: 6-tuple of fitted numbers
        isl_idx: island serial number
        g_idx: gaussian serial number
        flg: flagging (if any)
        """
        from . import functions as func
        import numpy as N

        # Add attribute definitions needed for output
        self.source_id_def = Int(doc="Source index", colname='Source_id')
        self.code_def = String(doc='Source code S, C, or M', colname='S_Code')
        self.gaus_num_def = Int(doc="Serial number of the gaussian for the image", colname='Gaus_id')
        self.island_id_def = Int(doc="Serial number of the island", colname='Isl_id')
        self.flag_def = Int(doc="Flag associated with gaussian", colname='Flag')
        self.total_flux_def = Float(doc="Total flux density, Jy", colname='Total_flux', units='Jy')
        self.total_fluxE_def = Float(doc="Total flux density error, Jy", colname='E_Total_flux',
                                     units='Jy')
        self.peak_flux_def = Float(doc="Peak flux density/beam, Jy/beam", colname='Peak_flux',
                                   units='Jy/beam')
        self.peak_fluxE_def = Float(doc="Peak flux density/beam error, Jy/beam",
                                    colname='E_Peak_flux', units='Jy/beam')
        self.centre_sky_def = List(Float(), doc="Sky coordinates of gaussian centre",
                                   colname=['RA', 'DEC'], units=['deg', 'deg'])
        self.centre_skyE_def = List(Float(), doc="Error on sky coordinates of gaussian centre",
                                    colname=['E_RA', 'E_DEC'], units=['deg', 'deg'])
        self.centre_pix_def = List(Float(), doc="Pixel coordinates of gaussian centre",
                                   colname=['Xposn', 'Yposn'], units=['pix', 'pix'])
        self.centre_pixE_def = List(Float(), doc="Error on pixel coordinates of gaussian centre",
                                    colname=['E_Xposn', 'E_Yposn'], units=['pix', 'pix'])
        self.size_sky_def = List(Float(), doc="Shape of the gaussian FWHM, PA, deg",
                                 colname=['Maj', 'Min', 'PA'], units=['deg', 'deg', 'deg'])
        self.size_skyE_def = List(Float(), doc="Error on shape of the gaussian FWHM, PA, deg",
                                  colname=['E_Maj', 'E_Min', 'E_PA'], units=['deg', 'deg', 'deg'])
        self.deconv_size_sky_def = List(Float(), doc="Deconvolved shape of the gaussian FWHM, PA, deg",
                                        colname=['DC_Maj', 'DC_Min', 'DC_PA'], units=['deg', 'deg', 'deg'])
        self.deconv_size_skyE_def = List(Float(), doc="Error on deconvolved shape of the gaussian FWHM, PA, deg",
                                         colname=['E_DC_Maj', 'E_DC_Min', 'E_DC_PA'], units=['deg', 'deg', 'deg'])
        self.size_sky_uncorr_def = List(Float(), doc="Shape in image plane of the gaussian FWHM, PA, deg",
                                        colname=['Maj_img_plane', 'Min_img_plane', 'PA_img_plane'],
                                        units=['deg', 'deg', 'deg'])
        self.size_skyE_uncorr_def = List(Float(), doc="Error on shape in image plane of the gaussian FWHM, PA, deg",
                                         colname=['E_Maj_img_plane', 'E_Min_img_plane', 'E_PA_img_plane'],
                                         units=['deg', 'deg', 'deg'])
        self.deconv_size_sky_uncorr_def = List(Float(), doc="Deconvolved shape in image plane of the gaussian FWHM, PA, deg",
                                               colname=['DC_Maj_img_plane', 'DC_Min_img_plane', 'DC_PA_img_plane'],
                                               units=['deg', 'deg', 'deg'])
        self.deconv_size_skyE_uncorr_def = List(Float(), doc="Error on deconvolved shape in image plane of the gaussian FWHM, PA, deg",
                                                colname=['E_DC_Maj_img_plane', 'E_DC_Min_img_plane', 'E_DC_PA_img_plane'],
                                                units=['deg', 'deg', 'deg'])
        self.rms_def = Float(doc="Island rms, Jy/beam", colname='Isl_rms', units='Jy/beam')
        self.mean_def = Float(doc="Island mean, Jy/beam", colname='Isl_mean', units='Jy/beam')
        self.total_flux_isl_def = Float(doc="Island total flux from sum of pixels", colname='Isl_Total_flux', units='Jy')
        self.total_flux_islE_def = Float(doc="Error on island total flux from sum of pixels", colname='E_Isl_Total_flux', units='Jy')
        self.gresid_rms_def = Float(doc="Island rms in Gaussian residual image", colname='Resid_Isl_rms', units='Jy/beam')
        self.gresid_mean_def = Float(doc="Island mean in Gaussian residual image", colname='Resid_Isl_mean', units='Jy/beam')
        self.sresid_rms_def = Float(doc="Island rms in Shapelet residual image", colname='Resid_Isl_rms', units='Jy/beam')
        self.sresid_mean_def = Float(doc="Island mean in Shapelet residual image", colname='Resid_Isl_mean', units='Jy/beam')
        self.wave_rms_def = Float(doc="Island rms in wavelet image, Jy/beam", colname='Wave_Isl_rms', units='Jy/beam')
        self.wave_mean_def = Float(doc="Island mean in wavelet image, Jy/beam", colname='Wave_Isl_mean', units='Jy/beam')
        self.jlevel_def = Int(doc="Wavelet number to which Gaussian belongs", colname='Wave_id')
        self.spec_indx_def = Float(doc="Spectral index", colname='Spec_Indx', units=None)
        self.e_spec_indx_def = Float(doc="Error in spectral index", colname='E_Spec_Indx', units=None)
        self.specin_flux_def = List(Float(), doc="Total flux density per channel, Jy", colname=['Total_flux'], units=['Jy'])
        self.specin_fluxE_def = List(Float(), doc="Error in total flux density per channel, Jy", colname=['E_Total_flux'], units=['Jy'])
        self.specin_freq_def = List(Float(), doc="Frequency per channel, Hz", colname=['Freq'], units=['Hz'])

        use_wcs = True
        self.gaussian_idx = g_idx
        self.gaus_num = 0  # stored later
        self.island_id = isl_idx
        self.jlevel = img.j
        self.flag = flg
        self.parameters = gaussian

        p = gaussian
        self.peak_flux = p[0]
        self.centre_pix = p[1:3]
        size = p[3:6]
        if func.approx_equal(size[0], img.pixel_beam()[0]*1.1) and \
                func.approx_equal(size[1], img.pixel_beam()[1]) and \
                func.approx_equal(size[2], img.pixel_beam()[2]+90.0) or \
                img.opts.fix_to_beam:
            # Check whether fitted Gaussian is just the distorted pixel beam given as an
            # initial guess (always set to [bm_maj*1.1, bm_min, bm_pa+90]) or if size was
            # fixed to the beam. If so, reset the size to the undistorted beam. Note:
            # these are sigma sizes, not FWHM sizes.
            size = img.pixel_beam()
            size = (size[0], size[1], size[2]+90.0)  # adjust angle so that corrected_size() works correctly
        size = func.corrected_size(size)  # gives fwhm and P.A.
        self.size_pix = size  # FWHM in pixels and P.A. CCW from +y axis

        # Use img.orig_beam for flux calculation and deconvolution on wavelet
        # images, as img.beam has been altered to match the wavelet scale.
        # Note: these are all FWHM sizes.
        if img.waveletimage:
            bm_pix = N.array(img.beam2pix(img.orig_beam))
        else:
            bm_pix = N.array(img.beam2pix(img.beam))

        # Calculate fluxes, sky sizes, etc. All sizes are FWHM.
        tot = p[0]*size[0]*size[1]/(bm_pix[0]*bm_pix[1])
        if flg == 0:
            # These are good Gaussians
            errors = func.get_errors(img, p+[tot], img.islands[isl_idx].rms, fixed_to_beam=img.opts.fix_to_beam)
            self.centre_sky = img.pix2sky(p[1:3])
            self.centre_skyE = img.pix2coord(errors[1:3], self.centre_pix, use_wcs=use_wcs)
            self.size_sky = img.pix2gaus(size, self.centre_pix, use_wcs=use_wcs)  # FWHM in degrees and P.A. east from north
            self.size_sky_uncorr = img.pix2gaus(size, self.centre_pix, use_wcs=False)  # FWHM in degrees and P.A. east from +y axis
            self.size_skyE = img.pix2gaus(errors[3:6], self.centre_pix, use_wcs=use_wcs, is_error=True)
            self.size_skyE_uncorr = img.pix2gaus(errors[3:6], self.centre_pix, use_wcs=False, is_error=True)
            gaus_dc, err = func.deconv2(bm_pix, size)
            self.deconv_size_sky = img.pix2gaus(gaus_dc, self.centre_pix, use_wcs=use_wcs)
            self.deconv_size_sky_uncorr = img.pix2gaus(gaus_dc, self.centre_pix, use_wcs=False)
            self.deconv_size_skyE = img.pix2gaus(errors[3:6], self.centre_pix, use_wcs=use_wcs, is_error=True)
            self.deconv_size_skyE_uncorr = img.pix2gaus(errors[3:6], self.centre_pix, use_wcs=False, is_error=True)
        else:
            # These are flagged Gaussians, so don't calculate sky values or errors
            errors = [0]*7
            self.centre_sky = [0., 0.]
            self.centre_skyE = [0., 0.]
            self.size_sky = [0., 0., 0.]
            self.size_sky_uncorr = [0., 0., 0.]
            self.size_skyE = [0., 0.]
            self.size_skyE_uncorr = [0., 0., 0.]
            self.deconv_size_sky = [0., 0., 0.]
            self.deconv_size_sky_uncorr = [0., 0., 0.]
            self.deconv_size_skyE = [0., 0., 0.]
            self.deconv_size_skyE_uncorr = [0., 0., 0.]
        self.total_flux = tot
        self.total_fluxE = errors[6]
        self.peak_fluxE = errors[0]
        self.total_fluxE = errors[6]
        self.centre_pixE = errors[1:3]
        self.size_pixE = errors[3:6]
        self.rms = img.islands[isl_idx].rms
        self.mean = img.islands[isl_idx].mean
        self.wave_rms = 0.0  # set if needed in the wavelet operation
        self.wave_mean = 0.0  # set if needed in the wavelet operation
        self.total_flux_isl = img.islands[isl_idx].total_flux
        self.total_flux_islE = img.islands[isl_idx].total_fluxE
