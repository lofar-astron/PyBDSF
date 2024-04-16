"""Module preprocess

Calculates some basic statistics of the image and sets up processing
parameters for PyBDSM.
"""
from __future__ import absolute_import

import numpy as N
from . import _cbdsm
from .image import *
from math import pi, sqrt, log
from . import const
from . import functions as func
from . import mylogger


class Op_preprocess(Op):
    """Preprocessing -- calculate some basic statistics and set
    processing parameters. Should assume that pixels outside the universe
    are blanked in QC ? """

    def __call__(self, img):
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Preprocess")
        bstat = func.bstat
        if img.opts.kappa_clip is None:
            kappa = -img.pixel_beamarea()
        else:
            kappa = img.opts.kappa_clip

        if img.opts.polarisation_do:
            pols = ['I', 'Q', 'U', 'V']
            ch0images = [img.ch0_arr, img.ch0_Q_arr, img.ch0_U_arr, img.ch0_V_arr]
            img.clipped_mean_QUV = []
            img.clipped_rms_QUV = []
        else:
            pols = ['I'] # assume I is always present
            ch0images = [img.ch0_arr]

        if hasattr(img, 'rms_mask'):
            mask = img.rms_mask
        else:
            mask = img.mask_arr
        opts = img.opts

        for ipol, pol in enumerate(pols):
            image = ch0images[ipol]

            ### basic stats
            mean, rms, cmean, crms, cnt = bstat(image, mask, kappa)
            if cnt > 198: cmean = mean; crms = rms
            if pol == 'I':
                if func.approx_equal(crms, 0.0, rel=None):
                    raise RuntimeError('Clipped rms appears to be zero. Check for regions '\
                                           'with values of 0 and\nblank them (with NaNs) '\
                                           'or use trim_box to exclude them.')
                img.raw_mean    = mean
                img.raw_rms     = rms
                img.clipped_mean= cmean
                img.clipped_rms = crms
                mylog.info('%s %.4f %s %.4f %s ' % ("Raw mean (Stokes I) = ", mean*1000.0, \
                           'mJy and raw rms = ',rms*1000.0, 'mJy'))
                mylog.info('%s %.4f %s %s %.4f %s ' % ("sigma clipped mean (Stokes I) = ", cmean*1000.0, \
                           'mJy and ','sigma clipped rms = ',crms*1000.0, 'mJy'))
            else:
                img.clipped_mean_QUV.append(cmean)
                img.clipped_rms_QUV.append(crms)
                mylog.info('%s %s %s %.4f %s %s %.4f %s ' % ("sigma clipped mean (Stokes ", pol, ") = ", cmean*1000.0, \
                           'mJy and ','sigma clipped rms = ',crms*1000.0, 'mJy'))

        image = img.ch0_arr
        # Check if pixels are outside the universe
        if opts.check_outsideuniv:
            mylogger.userinfo(mylog, "Checking for pixels outside the universe")
            noutside_univ = self.outside_univ(img)
            img.noutside_univ = noutside_univ
            frac_blank = round(float(noutside_univ)/float(image.shape[0]*image.shape[1]),3)
            mylogger.userinfo(mylog, "Number of additional pixels blanked", str(noutside_univ)
                              +' ('+str(frac_blank*100.0)+'%)')
        else:
            noutside_univ = 0

        # If needed, (re)mask the image
        if noutside_univ > 0:
            mask = N.isnan(img.ch0_arr)
            masked = mask.any()
            img.masked = masked
            if masked:
                img.mask_arr = mask
            img.blankpix = N.sum(mask)


        ### max/min pixel value & coordinates
        shape = image.shape[0:2]
        if mask is not None:
            img.blankpix = N.sum(mask)
        if img.blankpix == 0:
            max_idx = image.argmax()
            min_idx = image.argmin()
        else:
            max_idx = N.nanargmax(image)
            min_idx = N.nanargmin(image)

        img.maxpix_coord = N.unravel_index(max_idx, shape)
        img.minpix_coord = N.unravel_index(min_idx, shape)
        img.max_value    = image.flat[max_idx]
        img.min_value    = image.flat[min_idx]

        ### Solid angle of the image
        cdelt = N.array(img.wcs_obj.acdelt[:2])
        img.omega = N.product(shape)*abs(N.product(cdelt))/(180.*180./pi/pi)

        ### Total flux in ch0 image
        if 'atrous' in img.filename or img._pi or img.log == 'Detection image':
            # Don't do this estimate for atrous wavelet images
            # or polarized intensity image,
            # as it doesn't give the correct flux. Also, ignore
            # the flux in the detection image, as it's likely
            # wrong (e.g., not corrected for the primary beam).
            img.ch0_sum_jy = 0
        else:
            im_flux = N.nansum(image)/img.pixel_beamarea() # Jy
            img.ch0_sum_jy = im_flux
            mylogger.userinfo(mylog, 'Flux from sum of (non-blank) pixels',
                              '%.3f Jy' % (im_flux,))

        ### if image seems confused, then take background mean as zero instead
        alpha_sourcecounts = 2.5  # approx diff src count slope. 2.2?
        if opts.bmpersrc_th is None:
            if mask is not None:
                unmasked = N.where(~img.mask_arr)
                n = (image[unmasked] >= 5.*crms).sum()
            else:
                n = (image >= 5.*crms).sum()
            if n <= 0:
                n = 1
                mylog.info('No pixels in image > 5-sigma.')
                mylog.info('Taking number of pixels above 5-sigma as 1.')
            img.bmpersrc_th = N.product(shape)/((alpha_sourcecounts-1.)*n)
            mylog.info('%s %6.2f' % ('Estimated bmpersrc_th = ', img.bmpersrc_th))
        else:
            img.bmpersrc_th = opts.bmpersrc_th
            mylog.info('%s %6.2f' % ('Taking default bmpersrc_th = ', img.bmpersrc_th))

        confused = False
        if opts.mean_map == 'default':
            if img.bmpersrc_th <= 25. or cmean/crms >= 0.1:
                confused = True
        img.confused = confused
        mylog.info('Parameter confused is '+str(img.confused))

        img.completed_Ops.append('preprocess')
        return img

    def outside_univ(self,img):
        """ Checks if a pixel is outside the universe and is not blanked,
        and blanks it. (fits files written by CASA dont do this).  """

        noutside = 0
        n, m = img.ch0_arr.shape
        for i in range(n):
            for j in range(m):
                out = False
                err = ''
                pix1 = (i,j)
                try:
                    skyc = img.pix2sky(pix1)
                    pix2 = img.sky2pix(skyc)
                    if abs(pix1[0]-pix2[0]) > 0.5 or abs(pix1[1]-pix2[1]) > 0.5: out=True
                except RuntimeError as err:
                    pass
                if out or ("8" in str(err)):
                    noutside += 1
                    ch0 = img.ch0_arr
                    ch0[pix1] = float("NaN")
                    img.ch0_arr = ch0
        return noutside
