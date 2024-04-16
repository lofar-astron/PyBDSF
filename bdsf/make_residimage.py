"""Module make_residimage.

It calculates residual image from the list of gaussians and shapelets
"""
from __future__ import absolute_import

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

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"ResidImage")
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
        if isinstance(img.mask_arr, N.ndarray):
            pix_masked = N.where(img.mask_arr == True)
            model_gaus[pix_masked] = N.nan
            resid_gaus[pix_masked] = N.nan

        img.model_gaus_arr = model_gaus
        img.resid_gaus_arr = resid_gaus

        if img.opts.output_all or img.opts.savefits_residim:
            if img.waveletimage:
                resdir = img.basedir + '/wavelet/residual/'
            else:
                resdir = img.basedir + '/residual/'
            if not os.path.exists(resdir): os.makedirs(resdir)
            func.write_image_to_file(img.use_io, img.imagename + '.resid_gaus.fits', resid_gaus, img, resdir)
            mylog.info('%s %s' % ('Writing', resdir+img.imagename+'.resid_gaus.fits'))
        if img.opts.output_all or img.opts.savefits_modelim:
            if img.waveletimage:
                moddir = img.basedir + '/wavelet/model/'
            else:
                moddir = img.basedir + '/model/'
            if not os.path.exists(moddir): os.makedirs(moddir)
            func.write_image_to_file(img.use_io, img.imagename + '.model.fits', (img.ch0_arr - resid_gaus), img, moddir)
            mylog.info('%s %s' % ('Writing', moddir+img.imagename+'.model_gaus.fits'))

        ### residual rms and mean per island
        for isl in img.islands:
            resid = resid_gaus[tuple(isl.bbox)]
            self.calc_resid_mean_rms(isl, resid, type='gaus')

        # Calculate some statistics for the Gaussian residual image
        non_masked = N.where(~N.isnan(img.ch0_arr))
        mean = N.mean(resid_gaus[non_masked], axis=None)
        std_dev = N.std(resid_gaus[non_masked], axis=None)
        skew = stats.skew(resid_gaus[non_masked], axis=None)
        kurt = stats.kurtosis(resid_gaus[non_masked], axis=None)
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

            # Calculate some statistics for the Shapelet residual image
            non_masked = N.where(~N.isnan(img.ch0_arr))
            mean = N.mean(resid_shap[non_masked], axis=None)
            std_dev = N.std(resid_shap[non_masked], axis=None)
            skew = stats.skew(resid_shap[non_masked], axis=None)
            kurt = stats.kurtosis(resid_shap[non_masked], axis=None)
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
        """
        if len(isl.gaul) == 0:
            resid = N.zeros(isl.shape, dtype=N.float32)

        ind = N.where(~isl.mask_active)
        resid = resid[ind]
        if type == 'gaus':
            isl.gresid_rms = N.std(resid)
            isl.gresid_mean = N.mean(resid)
        else:
            isl.sresid_rms = N.std(resid)
            isl.sresid_mean = N.mean(resid)
        if hasattr(isl, 'sources'):
            for src in isl.sources:
                if type == 'gaus':
                    src.gresid_rms = N.std(resid)
                    src.gresid_mean = N.mean(resid)
                else:
                    src.sresid_rms = N.std(resid)
                    src.sresid_mean = N.mean(resid)
                for g in src.gaussians:
                    if type == 'gaus':
                        g.gresid_rms = N.std(resid)
                        g.gresid_mean = N.mean(resid)
                    else:
                        g.sresid_rms = N.std(resid)
                        g.sresid_mean = N.mean(resid)
        if hasattr(isl, 'dsources'):
            for dsrc in isl.dsources: # Handle dummy sources (if any)
                if type == 'gaus':
                    dsrc.gresid_rms = N.std(resid)
                    dsrc.gresid_mean = N.mean(resid)
                else:
                    dsrc.sresid_rms = N.std(resid)
                    dsrc.sresid_mean = N.mean(resid)
