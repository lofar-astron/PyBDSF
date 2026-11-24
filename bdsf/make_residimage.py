"""Module make_residimage.

It calculates residual image from the list of gaussians and shapelets
"""

import numpy as N
from scipy import stats # for skew and kurtosis
from .image import *
from .shapelets import *
from . import mylogger


class Op_make_residimage(Op):
    """Creates an image from the fitted gaussians
    or shapelets.

    The resulting model image is stored in the
    resid_gaus or resid_shap attribute.

    Prerequisites: module gausfit or shapelets should
    be run first.
    """

    def __call__(self, img):
        from . import functions as func
        from copy import deepcopy as cp
        import os

        mylog = mylogger.logging.getLogger("PyBDSF."+img.log+"ResidImage")
        mylog.info("Calculating residual image after subtracting reconstructed gaussians")
        shape = img.ch0_arr.shape
        thresh= img.opts.fittedimage_clip

        resid_gaus = cp(img.ch0_arr)
        model_gaus = N.zeros(shape, dtype=N.float32)
        for g in img.gaussians:
            C1, C2 = g.centre_pix
            if hasattr(g, 'wisland_id') and img.waveletimage:
                isl = img.islands[g.wisland_id]
            else:
                isl = img.islands[g.island_id]
            b = self.find_bbox(thresh*isl.rms, g)

            bbox = N.s_[max(0, int(C1-b)):min(shape[0], int(C1+b+1)),
                        max(0, int(C2-b)):min(shape[1], int(C2+b+1))]

            x_ax, y_ax = N.mgrid[bbox]
            ffimg = func.gaussian_fcn(g, x_ax, y_ax)
            resid_gaus[bbox] = resid_gaus[bbox] - ffimg
            model_gaus[bbox] = model_gaus[bbox] + ffimg

        # Apply mask to model and resid images
        if hasattr(img, 'rms_mask'):
            mask = img.rms_mask
        else:
            mask = img.mask_arr
        if isinstance(mask, N.ndarray):
            pix_masked = N.where(mask == True)
            model_gaus[pix_masked] = N.nan
            resid_gaus[pix_masked] = N.nan

        img.model_gaus_arr = model_gaus
        img.resid_gaus_arr = resid_gaus

        if img.opts.resid_sigma_stats_do and not img.waveletimage and not img._pi and \
            not (img.opts.atrous_do and 'wavelet_atrous' in img.completed_Ops):
            # Use the same valid pixels for before and after statistics.
            r_b_all = N.where(N.isfinite(resid_gaus), img.ch0_arr, N.nan)
            s1_b_all, ss_b_all, n_all   = func.residual_sumsq_in_sigma(r_b_all,    img.mean_arr, img.rms_arr)
            s1_a_all, ss_a_all, n_all_a = func.residual_sumsq_in_sigma(resid_gaus, img.mean_arr, img.rms_arr)
            rms_before = N.sqrt(ss_b_all / n_all) if n_all > 0 else N.nan
            rms_after = N.sqrt(ss_a_all / n_all_a) if n_all_a > 0 else N.nan
            img._resid_gaus_rms_sigma = rms_after
            mylogger.userinfo(mylog, 'Calculating sigma statistics of residual image (\'resid_sigma_stats_do\')')
            mylogger.userinfo(mylog, "Residual RMS  before / after Gaussian fitting ", "%.3f / %.3f σ" % (rms_before, rms_after))
            if N.isfinite(rms_before) and rms_before > 0 and N.isfinite(rms_after):
                change = 100.0 * (rms_after / rms_before - 1.0)
                mylogger.userinfo(mylog, "Residual RMS  change after Gaussian fitting   ", "%+.2f%%" % change)

            # The same statistic separately inside and outside islands
            s1_b_in = s1_a_in = 0.0
            ss_b_in = ss_a_in = 0.0
            n_in = 0
            for isl in img.islands:
                bb = tuple(isl.bbox)
                # Active island pixels only; NaN elsewhere is skipped by the helper.
                # 'before' uses the same pixel set as 'after' (also excludes masked pixels).
                r_a = N.where(~isl.mask_active, resid_gaus[bb], N.nan)
                r_b = N.where(N.isfinite(r_a), img.ch0_arr[bb], N.nan)
                m = img.mean_arr[bb]
                s = img.rms_arr[bb]
                s1b, sb, _  = func.residual_sumsq_in_sigma(r_b, m, s)
                s1a, sa, na = func.residual_sumsq_in_sigma(r_a, m, s)
                s1_b_in += s1b
                s1_a_in += s1a
                ss_b_in += sb
                ss_a_in += sa
                n_in += na
            n_out = n_all_a - n_in
            if n_in > 0 and n_out > 0:
                mylogger.userinfo(mylog, "Residual RMS  inside  islands (before / after)", "%.3f / %.3f σ"   % (N.sqrt(ss_b_in / n_in),               N.sqrt(ss_a_in / n_in)))
                mylogger.userinfo(mylog, "Residual RMS  outside islands (before / after)", "%.3f / %.3f σ"   % (N.sqrt((ss_b_all - ss_b_in) / n_out), N.sqrt((ss_a_all - ss_a_in) / n_out)))
                mylogger.userinfo(mylog, "Residual mean inside  islands (before / after)", "%+.3f / %+.3f σ" % (s1_b_in / n_in,                       s1_a_in / n_in))
                mylogger.userinfo(mylog, "Residual mean outside islands (before / after)", "%+.3f / %+.3f σ" % ((s1_b_all - s1_b_in) / n_out,         (s1_a_all - s1_a_in) / n_out))
            if not img.opts.residual_stats_do:
                mylog.info("More residual statistics can be enabled with the 'residual_stats_do' option.")

        if img.opts.output_all or img.opts.savefits_residim:
            if img.waveletimage:
                resdir = img.basedir + '/wavelet/residual/'
            else:
                resdir = img.basedir + '/residual/'
            os.makedirs(resdir, exist_ok=True)
            func.write_image_to_file(img.use_io, img.imagename + '.resid_gaus.fits', resid_gaus, img, resdir)
            mylog.info('%s %s' % ('Writing', resdir+img.imagename+'.resid_gaus.fits'))
        if img.opts.output_all or img.opts.savefits_modelim:
            if img.waveletimage:
                moddir = img.basedir + '/wavelet/model/'
            else:
                moddir = img.basedir + '/model/'
            os.makedirs(moddir, exist_ok=True)
            func.write_image_to_file(img.use_io, img.imagename + '.model_gaus.fits', (img.ch0_arr - resid_gaus), img, moddir)
            mylog.info('%s %s' % ('Writing', moddir+img.imagename+'.model_gaus.fits'))

        ### residual rms and mean per island
        no_gaus_ids = []
        for isl in img.islands:
            if len(isl.gaul) == 0:
                no_gaus_ids.append(isl.island_id)
            resid = resid_gaus[tuple(isl.bbox)]
            self.calc_resid_mean_rms(isl, resid, type='gaus')

        if no_gaus_ids and img.opts.verbose_fitting:
            shown = ', '.join(str(i) for i in no_gaus_ids[:20])
            if len(no_gaus_ids) > 20:
                shown += ', ...'
            mylog.warning("%d island(s) have no fitted Gaussians; their residual "
                          "mean/rms are computed from the unmodified image data "
                          "(island ids: %s)" % (len(no_gaus_ids), shown))

        if img.opts.residual_stats_do:
            # Calculate some statistics for the Gaussian residual image
            resid_gaus_non_masked = resid_gaus[~N.isnan(resid_gaus)]
            mean = N.mean(resid_gaus_non_masked, axis=None)
            std_dev = N.std(resid_gaus_non_masked, axis=None)
            skew = stats.skew(resid_gaus_non_masked, axis=None)
            kurt = stats.kurtosis(resid_gaus_non_masked, axis=None)
            stat_msg = "Statistics of the Gaussian residual image:\n"
            stat_msg += "        mean: %.3e (Jy/beam)\n" % mean
            stat_msg += "    std. dev: %.3e (Jy/beam)\n" % std_dev
            stat_msg += "        skew: %.3f\n" % skew
            stat_msg += "    kurtosis: %.3f" % kurt
            mylog.info(stat_msg)

        # Now residual image for shapelets
        if img.opts.shapelet_do:
            mylog.info("Calculating residual image after subtracting reconstructed shapelets")
            shape = img.ch0_arr.shape
            fimg = N.zeros(shape, dtype=N.float32)

            for isl in img.islands:
                if  hasattr(isl, 'shapelet_beta'):
                    if isl.shapelet_beta > 0: # make sure shapelet has nonzero scale for this island
                        mask=isl.mask_active
                        cen=isl.shapelet_centre-N.array(isl.origin)
                        basis, beta, nmax, cf = isl.shapelet_basis, isl.shapelet_beta, \
                                                isl.shapelet_nmax, isl.shapelet_cf
                        image_recons=reconstruct_shapelets(isl.shape, mask, basis, beta, cen, nmax, cf)
                        fimg[tuple(isl.bbox)] += image_recons

            model_shap = fimg
            resid_shap = img.ch0_arr - fimg
            if img.opts.shapelet_gresid:
                # also subtract Gaussian model image
                shape = img.ch0_arr.shape
                thresh= img.opts.fittedimage_clip
                model_gaus = N.zeros(shape, dtype=N.float32)
                for isl in img.islands:
                    for g in isl.gaul:
                        C1, C2 = g.centre_pix
                        b = self.find_bbox(thresh*isl.rms, g)
                        bbox = N.s_[max(0, int(C1-b)):min(shape[0], int(C1+b+1)),
                                    max(0, int(C2-b)):min(shape[1], int(C2+b+1))]
                        x_ax, y_ax = N.mgrid[bbox]
                        ffimg = func.gaussian_fcn(g, x_ax, y_ax)
                        model_gaus[bbox] = model_gaus[bbox] + ffimg
                resid_shap -= model_gaus

            # Apply mask to model and resid images
            if hasattr(img, 'rms_mask'):
                mask = img.rms_mask
            else:
                mask = img.mask_arr
            if isinstance(mask, N.ndarray):
                pix_masked = N.where(mask == True)
                model_shap[pix_masked] = N.nan
                resid_shap[pix_masked] = N.nan

            img.model_shap_arr = model_shap
            img.resid_shap_arr = resid_shap

            if img.opts.output_all:
                func.write_image_to_file(img.use_io, img.imagename + '.resid_shap.fits', resid_shap, img, resdir)
                mylog.info('%s %s' % ('Writing ', resdir+img.imagename+'.resid_shap.fits'))

            ### shapelet residual rms and mean per island
            for isl in img.islands:
                resid = resid_shap[tuple(isl.bbox)]
                self.calc_resid_mean_rms(isl, resid, type='shap')

            if img.opts.residual_stats_do:
                # Calculate some statistics for the Shapelet residual image
                resid_shap_non_masked = resid_shap[~N.isnan(resid_shap)]
                mean = N.mean(resid_shap_non_masked, axis=None)
                std_dev = N.std(resid_shap_non_masked, axis=None)
                skew = stats.skew(resid_shap_non_masked, axis=None)
                kurt = stats.kurtosis(resid_shap_non_masked, axis=None)
                mylog.info("Statistics of the Shapelet residual image:")
                mylog.info("        mean: %.3e (Jy/beam)" % mean)
                mylog.info("    std. dev: %.3e (Jy/beam)" % std_dev)
                mylog.info("        skew: %.3f" % skew)
                mylog.info("    kurtosis: %.3f" % kurt)

        img.completed_Ops.append('make_residimage')
        return img

    def find_bbox(self, thresh, g):
        """Calculate bounding box for gaussian.

        This function calculates size of the box for evaluating
        gaussian, so that value of gaussian is smaller than threshold
        outside of the box.

        Parameters:
        thres: threshold
        g: Gaussian object
        """

        from math import ceil, sqrt, log
        A = g.peak_flux
        S = g.size_pix[0]
        if A == 0.0:
            return ceil(S*1.5)
        if thresh/A >= 1.0 or thresh/A <= 0.0:
            return ceil(S*1.5)
        return ceil(S*sqrt(-2*log(thresh/A)))

    def calc_resid_mean_rms(self, isl, resid, type):
        """Inserts mean and rms of residual image into isl, src, and gaussians

        type - specifies 'gaus' or 'shap'

        If an island has no fitted Gaussians, the residual data are left
        unaltered (nothing was subtracted), so the statistics describe the
        local background noise plus any unmodeled emission.
        """

        ind = N.where(~isl.mask_active)
        resid = resid[ind]
        if type == 'gaus':
            isl.gresid_rms = N.nanstd(resid)
            isl.gresid_mean = N.nanmean(resid)
        else:
            isl.sresid_rms = N.nanstd(resid)
            isl.sresid_mean = N.nanmean(resid)
        if hasattr(isl, 'sources'):
            for src in isl.sources:
                if type == 'gaus':
                    src.gresid_rms = N.nanstd(resid)
                    src.gresid_mean = N.nanmean(resid)
                else:
                    src.sresid_rms = N.nanstd(resid)
                    src.sresid_mean = N.nanmean(resid)
                for g in src.gaussians:
                    if type == 'gaus':
                        g.gresid_rms = N.nanstd(resid)
                        g.gresid_mean = N.nanmean(resid)
                    else:
                        g.sresid_rms = N.nanstd(resid)
                        g.sresid_mean = N.nanmean(resid)
        if hasattr(isl, 'dsources'):
            for dsrc in isl.dsources: # Handle dummy sources (if any)
                if type == 'gaus':
                    dsrc.gresid_rms = N.nanstd(resid)
                    dsrc.gresid_mean = N.nanmean(resid)
                else:
                    dsrc.sresid_rms = N.nanstd(resid)
                    dsrc.sresid_mean = N.nanmean(resid)
