"""Module make_residimage.

It calculates residual image from the list of gaussians and shapelets
"""

import numpy as N
from scipy import stats # for skew and kurtosis
from image import *
from shapelets import *
import mylogger

### Insert attribute into Image class for model image
Image.resid_gaus = NArray(doc="Residual image calculated from " \
                                "extracted gaussians")
Image.resid_shap = NArray(doc="Residual image calculated from " \
                                "shapelet coefficient")
Image.model_gaus = NArray(doc="Model image calculated from " \
                                "extracted gaussians")
Image.model_shap = NArray(doc="Model image calculated from " \
                                "shapelet coefficient")

class Op_make_residimage(Op):
    """Creates an image from the fitted gaussians
    or shapelets.

    The resulting model image is stored in the
    resid_gaus or resid_shap attribute.

    Prerequisites: module gausfit or shapelets should
    be run first.
    """

    def __call__(self, img):
        import functions as func
        from copy import deepcopy as cp
        import os

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"ResidImage")
        mylog.info("Calculating residual image after subtracting reconstructed gaussians")
        shape = img.ch0.shape
        thresh= img.opts.fittedimage_clip

        img.resid_gaus = cp(img.ch0)
        img.model_gaus = N.zeros(shape, dtype=float)
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
            img.resid_gaus[bbox] = img.resid_gaus[bbox] - ffimg
            img.model_gaus[bbox] = img.model_gaus[bbox] + ffimg

        # Apply mask to model and resid images
        if hasattr(img, 'rms_mask'):
            mask = img.rms_mask
        else:
            mask = img.mask
        if isinstance(img.mask, N.ndarray):
            pix_masked = N.where(img.mask == True)
            img.model_gaus[pix_masked] = N.nan
            img.resid_gaus[pix_masked] = N.nan

        if img.opts.output_all:
            if img.waveletimage:
                resdir = img.basedir + '/wavelet/residual/'
                moddir = img.basedir + '/wavelet/model/'
            else:
                resdir = img.basedir + '/residual/'
                moddir = img.basedir + '/model/'
            if not os.path.exists(resdir): os.mkdir(resdir)
            if not os.path.exists(moddir): os.mkdir(moddir)
            func.write_image_to_file(img.use_io, img.imagename + '.resid_gaus.fits', img.resid_gaus, img, resdir)
            mylog.info('%s %s' % ('Writing', resdir+img.imagename+'.resid_gaus.fits'))
            func.write_image_to_file(img.use_io, img.imagename + '.model.fits', (img.ch0 - img.resid_gaus), img, moddir)
            mylog.info('%s %s' % ('Writing', moddir+img.imagename+'.model_gaus.fits'))

        ### residual rms and mean per island
        for isl in img.islands:
            resid = img.resid_gaus[isl.bbox]
            n, m = resid.shape

            ind = N.where(~isl.mask_active)
            resid = resid[ind]
            isl.gresid_rms = N.std(resid)
            isl.gresid_mean = N.mean(resid)
            for src in isl.sources:
                src.gresid_rms = N.std(resid)
                src.gresid_mean = N.mean(resid)
                for g in src.gaussians:
                    g.gresid_rms = N.std(resid)
                    g.gresid_mean = N.mean(resid)

        # Calculate some statistics for the Gaussian residual image
        non_masked = N.where(~N.isnan(img.ch0))
        mean = N.mean(img.resid_gaus[non_masked], axis=None)
        std_dev = N.std(img.resid_gaus[non_masked], axis=None)
        skew = stats.skew(img.resid_gaus[non_masked], axis=None)
        kurt = stats.kurtosis(img.resid_gaus[non_masked], axis=None)
        mylog.info("Statistics of the Gaussian residual image:")
        mylog.info("        mean: %.3e (Jy/beam)" % mean)
        mylog.info("    std. dev: %.3e (Jy/beam)" % std_dev)
        mylog.info("        skew: %.3f" % skew)
        mylog.info("    kurtosis: %.3f" % kurt)

        # Now residual image for shapelets
        if img.opts.shapelet_do:
            mylog.info("Calculating residual image after subtracting reconstructed shapelets")
            shape = img.ch0.shape
            fimg = N.zeros(shape, dtype=float)

            for isl in img.islands:
              if isl.shapelet_beta > 0: # make sure shapelet has nonzero scale for this island
                mask=isl.mask_active
                cen=isl.shapelet_centre-N.array(isl.origin)
                basis, beta, nmax, cf = isl.shapelet_basis, isl.shapelet_beta, \
                                        isl.shapelet_nmax, isl.shapelet_cf
                image_recons=reconstruct_shapelets(isl.shape, mask, basis, beta, cen, nmax, cf)
                fimg[isl.bbox] += image_recons

            img.model_shap = fimg
            img.resid_shap = img.ch0 - fimg
            # Apply mask to model and resid images
            if hasattr(img, 'rms_mask'):
                mask = img.rms_mask
            else:
                mask = img.mask
            if isinstance(mask, N.ndarray):
                pix_masked = N.where(mask == True)
                img.model_shap[pix_masked] = N.nan
                img.resid_shap[pix_masked] = N.nan

            if img.opts.output_all:
                func.write_image_to_file(img.use_io, img.imagename + '.resid_shap.fits', img.resid_shap, img, resdir)
                mylog.info('%s %s' % ('Writing ', resdir+img.imagename+'.resid_shap.fits'))

            ### shapelet residual rms and mean per island
            for isl in img.islands:
                resid = img.resid_shap[isl.bbox]
                n, m = resid.shape
                ind = N.where(~isl.mask_active)
                resid = resid[ind]
                isl.sresid_rms = N.std(resid)
                isl.sresid_mean = N.mean(resid)
                for src in isl.sources:
                    src.sresid_rms = N.std(resid)
                    src.sresid_mean = N.mean(resid)
                    for g in src.gaussians:
                        g.sresid_rms = N.std(resid)
                        g.sresid_mean = N.mean(resid)

            # Calculate some statistics for the Shapelet residual image
            non_masked = N.where(~N.isnan(img.ch0))
            mean = N.mean(img.resid_gaus[non_masked], axis=None)
            std_dev = N.std(img.resid_gaus[non_masked], axis=None)
            skew = stats.skew(img.resid_gaus[non_masked], axis=None)
            kurt = stats.kurtosis(img.resid_gaus[non_masked], axis=None)
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


