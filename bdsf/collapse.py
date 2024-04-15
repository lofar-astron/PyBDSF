"""Module collapse

Defines operation Op_collapse which collapses 3D image. Calculates and
stores mean and rms (normal and clipped) per channel anyway for further
use, even if weights are unity.
"""
from __future__ import absolute_import

import numpy as N
from .image import *
from . import _cbdsm
#_cbdsm.init_numpy()
from . import mylogger
from . import functions as func


class Op_collapse(Op):
    """Collapse 3D image"""

    def __call__(self, img):
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Collapse")
        if img.opts.polarisation_do:
            pols = ['I', 'Q', 'U', 'V'] # make sure I is done first
        else:
            pols = ['I'] # assume I is always present
            img.ch0_Q_arr = None
            img.ch0_U_arr = None
            img.ch0_V_arr = None

        if img.shape[1] > 1:
            c_mode = img.opts.collapse_mode
            chan0 = img.opts.collapse_ch0
            c_list = img.opts.collapse_av
            c_wts = img.opts.collapse_wt
            if not c_list:
                c_list = N.arange(img.shape[1])
            if len(c_list) == 1 and c_mode=='average':
                c_mode = 'single'
                chan0 = c_list[0]
                img.collapse_ch0 = chan0
            ch0sh = img.image_arr.shape[2:]
            if img.opts.polarisation_do:
                ch0images = ['ch0_arr', 'ch0_Q_arr', 'ch0_U_arr', 'ch0_V_arr']
            else:
                ch0images = ['ch0_arr']

            # Check whether the collapse channel index is sensible
            if chan0 < 0 or chan0 >= len(c_list):
                raise RuntimeError('The channel index (set with the "collapse_ch0" option) '
                                   'must be greater than zero and less than the number of '
                                   'channels ({}).'.format(len(c_list)))

            # assume all Stokes images have the same blank pixels as I:
            blank = N.isnan(img.image_arr[0])
            hasblanks = blank.any()
            if img.opts.kappa_clip is None:
                kappa = -img.pixel_beamarea()
            else:
                kappa = img.opts.kappa_clip

            mean, rms, cmean, crms = chan_stats(img, kappa)
            img.channel_mean = mean; img.channel_rms = rms
            img.channel_clippedmean = cmean; img.channel_clippedrms = crms

            for ipol, pol in enumerate(pols):
                if c_mode == 'single':
                    if pol == 'I':
                        ch0 = img.image_arr[0, chan0]
                        img.ch0_arr = ch0

                        # Construct weights so that desired channel has weight of 1 and all
                        # others have weight of 0.  The init_freq_collapse function will then
                        # select the intended frequency
                        wtarr = N.zeros(len(c_list))
                        wtarr[chan0] = 1.0
                        init_freq_collapse(img, wtarr)
                        mylogger.userinfo(mylog, 'Source extraction will be ' \
                                              'done on channel', '%i (%.3f MHz)' % \
                                              (chan0, img.frequency/1e6))
                    else:
                        ch0[:] = img.image_arr[ipol, chan0][:]
                        img.__setattr__(ch0images[ipol][:], ch0)

                elif c_mode == 'average':
                    if not hasblanks:
                        if pol == 'I':
                            ch0, wtarr = avspc_direct(c_list, img.image_arr[0], img.channel_clippedrms, c_wts)
                        else:
                            # use wtarr from the I image, which is always collapsed first
                            ch0, wtarr = avspc_direct(c_list, img.image_arr[ipol], img.channel_clippedrms, c_wts, wtarr=wtarr)
                    else:
                        if pol == 'I':
                            ch0, wtarr = avspc_blanks(c_list, img.image_arr[0], img.channel_clippedrms, c_wts)
                        else:
                            # use wtarr from the I image, which is always collapsed first
                            ch0, wtarr = avspc_blanks(c_list, img.image_arr[ipol], img.channel_clippedrms, c_wts, wtarr=wtarr)
                    img.__setattr__(ch0images[ipol][:], ch0)

                    if pol == 'I':
                        img.avspc_wtarr = wtarr
                        init_freq_collapse(img, wtarr)
                        if c_wts == 'unity':
                            mylogger.userinfo(mylog, 'Channels averaged with '\
                                                  'uniform weights')
                        else:
                            mylogger.userinfo(mylog, 'Channels averaged with '\
                                                  'weights=(1/rms)^2')
                        mylogger.userinfo(mylog, 'Source extraction will be '\
                                              'done on averaged ("ch0") image')
                        mylogger.userinfo(mylog, 'Frequency of averaged '\
                                              'image', '%.3f MHz' % \
                                              (img.frequency/1e6,))
                        str1 = " ".join(str(n) for n in c_list)
                        mylog.debug('%s %s' % ('Channels averaged : ', str1))
                        str1 = " ".join(["%9.4e" % n for n in wtarr])
                        mylog.debug('%s %s %s' % ('Channel weights : ', str1, '; unity=zero if c_wts="rms"'))
                elif c_mode=='file':
                    mylogger.userinfo(mylog, 'Reading ch0 image from file %s' % (img.opts.collapse_file))
                    image,hdr=func.read_image_from_file(img.opts.collapse_file, img, None, quiet=False)
                    if pol == 'I':
                        ch0 = image[0,0]
                        img.ch0_arr = ch0

                    else:
                        raise NotImplementedError('Polarization cubes not allowed in file mode')
                else:
                    raise NotImplementedError('Mode supplied not implemented') # should never happen!
                if img.opts.output_all:
                    func.write_image_to_file(img.use_io, img.imagename+'.ch0_'+pol+'.fits', ch0,
                                             img, outdir=img.basedir)
                    mylog.debug('%s %s ' % ('Writing file ', img.imagename+'.ch0_'+pol+'.fits'))

        else:
                # Only one channel in image
            image = img.image_arr
            img.ch0_arr = image[0, 0]
            mylogger.userinfo(mylog, 'Frequency of image',
                              '%.3f MHz' % (img.frequency/1e6,))
            if img.opts.polarisation_do:
                for pol in pols[1:]:
                    if pol == 'Q':
                        img.ch0_Q_arr = image[1, 0][:]
                    if pol == 'U':
                        img.ch0_U_arr = image[2, 0][:]
                    if pol == 'V':
                        img.ch0_V_arr = image[3, 0][:]

        # create mask if needed (assume all pols have the same mask as I)
        image = img.ch0_arr
        mask = N.isnan(image)
        img.blankpix = N.sum(mask)
        frac_blank = round(
            float(img.blankpix) / float(image.shape[0] * image.shape[1]),
            3)
        mylogger.userinfo(mylog, "Number of blank pixels", str(img.blankpix)
                          + ' (' + str(frac_blank * 100.0) + '%)')

        if img.opts.blank_limit is not None:
            import scipy
            import sys
            threshold = img.opts.blank_limit
            mylogger.userinfo(mylog, "Blanking pixels with values "
                              "below %.1e Jy/beam" % (threshold,))
            bad = (abs(image) < threshold)
            original_stdout = sys.stdout  # keep a reference to STDOUT
            sys.stdout = func.NullDevice()  # redirect the real STDOUT
            count = scipy.signal.convolve2d(bad, N.ones((3, 3)), mode='same')
            sys.stdout = original_stdout  # turn STDOUT back on
            mask_low = (count >= 5)
            image[N.where(mask_low)] = N.nan
            mask = N.isnan(image)
            img.blankpix = N.sum(mask)
            frac_blank = round(
                float(img.blankpix) / float(image.shape[0] *
                image.shape[1]), 3)
            mylogger.userinfo(mylog, "Total number of blanked pixels",
                str(img.blankpix) + ' (' + str(frac_blank * 100.0) + '%)')

        masked = mask.any()
        img.masked = masked
        if masked:
            img.mask_arr = mask
        else:
            img.mask_arr = None

        if img.blankpix == image.shape[0] * image.shape[1]:
            # ALL pixels are blanked!
            raise RuntimeError('All pixels in the image are blanked.')

        img.completed_Ops.append('collapse')


########################################################################################

def chan_stats(img, kappa):

    bstat = func.bstat #_cbdsm.bstat
    nchan = img.shape[1]
    mean = []; rms = []; cmean = []; crms = []
    for ichan in range(nchan):
        if isinstance(img, Image): # check if img is an Image or just an ndarray
            im = img.image_arr[0, ichan]
        else:
            im = img[0, ichan]

        if N.any(im):
            immask = N.isnan(im)
            if immask.all():
                m, r, cm, cr = 0, 0, 0, 0
            else:
                if immask.any():
                    m, r, cm, cr, cnt = bstat(im, immask, kappa)
                else:
                    m, r, cm, cr, cnt = bstat(im, None, kappa)
        else:
            m, r, cm, cr = 0, 0, 0, 0
        mean.append(m); rms.append(r); cmean.append(cm); crms.append(cr)

    return N.array(mean), N.array(rms), N.array(cmean), N.array(crms)


########################################################################################

def avspc_direct(c_list, image, rmsarr, c_wts, wtarr=None):

    shape2 = image.shape[1:]
    ch0 = N.zeros(shape2, dtype=N.float32)
    sumwts = 0.0
    if wtarr is None:
        wtarr = N.zeros(len(c_list))
        for i, ch in enumerate(c_list):
            im = image[ch]
            r = rmsarr[ch]
            if c_wts == 'unity': wt = 1.0
            if c_wts == 'rms': wt = r
            if r != 0:
                wt = 1.0/(wt*wt)
            else:
                wt = 0
            sumwts += wt
            ch0 += im*wt
            wtarr[i] = wt
    else:
        for i, ch in enumerate(c_list):
            im = image[ch]
            sumwts += wtarr[i]
            ch0 += im*wtarr[i]
    ch0 = ch0/sumwts

    return ch0, wtarr

########################################################################################

def avspc_blanks(c_list, image, rmsarr, c_wts, wtarr=None):

    shape2 = image.shape[1:]
    ch0 = N.zeros(shape2, dtype=N.float32)
    sumwtim = N.zeros(shape2, dtype=N.float32)
    if wtarr is None:
        wtarr = N.zeros(len(c_list))
        for i, ch in enumerate(c_list):
            im = image[ch]
            r = rmsarr[ch]
            if c_wts == 'unity': wt = 1.0
            if c_wts == 'rms': wt = r
            if r > 1e-18 and r < 1e18:
                # Set reasonable limits to avoid overflow of float32
                wt = 1.0/(wt*wt)
            else:
                wt = 0
            wtim = N.ones(shape2, dtype=N.float32)*wt*(~N.isnan(im))
            sumwtim += wtim
            ch0 += N.nan_to_num(im)*wtim
            wtarr[i] = wt
    else:
        for i, ch in enumerate(c_list):
            im = image[ch]
            wtim = N.ones(shape2)*wtarr[i]*(~N.isnan(im))
            sumwtim += wtim
            ch0 += N.nan_to_num(im)*wtim
    ch0 = ch0/sumwtim

    return ch0, wtarr

########################################################################################

def init_freq_collapse(img, wtarr):
    # Place appropriate, post-collapse frequency info in img
    # Calculate weighted average frequency
    if img.opts.frequency_sp is not None:
        c_list = img.opts.collapse_av
        if c_list == []: c_list = N.arange(img.image_arr.shape[1])
        freqs = img.opts.frequency_sp
        if len(freqs) != len(c_list):
            raise RuntimeError("Number of channels and number of frequencies specified "\
                         "by user do not match")
        sumwts = 0.0
        sumfrq = 0.0
        for i, ch in enumerate(c_list):
            sumwts += wtarr[i]
            sumfrq += freqs[ch]*wtarr[i]
        img.frequency = sumfrq / sumwts
        img.freq_pars = (img.frequency, 0.0, 0.0)
    else:
        # Calculate from header info
        c_list = img.opts.collapse_av
        if c_list == []: c_list = N.arange(img.image_arr.shape[1])
        sumwts = 0.0
        sumfrq = 0.0
        spec_indx = img.wcs_obj.wcs.spec
        if spec_indx == -1 and img.opts.frequency_sp is None:
            raise RuntimeError("Frequency information not found in header and frequencies "\
                         "not specified by user")
        else:
            for i, ch in enumerate(c_list):
                sumwts += wtarr[i]
                freq = img.wcs_obj.p2f(ch)
                sumfrq += freq*wtarr[i]
            img.frequency = sumfrq / sumwts
