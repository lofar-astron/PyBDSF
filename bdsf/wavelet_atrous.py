"""Compute a-trous wavelet transform of the gaussian residual image.

Do source extraction on this if asked.
"""
from __future__ import print_function
from __future__ import absolute_import
import numpy as N
from .image import *
from . import mylogger
import os
from . import has_pl
if has_pl:
    import matplotlib.pyplot as pl
from math import log, floor, sqrt
from .const import fwsig
from copy import deepcopy as cp
from . import functions as func
import gc
from numpy import array, product
import scipy.signal
from .preprocess import Op_preprocess
from .rmsimage import Op_rmsimage
from .threshold import Op_threshold
from .islands import Op_islands
from .gausfit import Op_gausfit, Gaussian
from .gaul2srl import Op_gaul2srl
from .make_residimage import Op_make_residimage
from .interface import raw_input_no_history
from . import statusbar
try:
    import pyfftw.interfaces
    pyfftw.interfaces.cache.enable()
    N.fft.fftn = pyfftw.interfaces.numpy_fft.fftn
    N.fft.ifftn = pyfftw.interfaces.numpy_fft.ifftn
    scipy.signal.signaltools.fftn = pyfftw.interfaces.scipy_fftpack.fftn
    scipy.signal.signaltools.ifftn = pyfftw.interfaces.scipy_fftpack.ifftn
    has_pyfftw = True
except ImportError:
    has_pyfftw = False


class Op_wavelet_atrous(Op):
    """Compute a-trous wavelet transform of the gaussian residual image."""

    def __call__(self, img):

        mylog = mylogger.logging.getLogger("PyBDSF." + img.log + "Wavelet")

        if img.opts.atrous_do:
            if img.nisl == 0:
                mylog.warning("No islands found. Skipping wavelet decomposition.")
                img.completed_Ops.append('wavelet_atrous')
                return

            mylog.info("Decomposing gaussian residual image into a-trous wavelets")
            bdir = img.basedir + '/wavelet/'
            if img.opts.output_all:
                if not os.path.isdir(bdir):
                    os.makedirs(bdir)
                if not os.path.isdir(bdir + '/residual/'):
                    os.makedirs(bdir + '/residual/')
                if not os.path.isdir(bdir + '/model/'):
                    os.makedirs(bdir + '/model/')
            dobdsm = img.opts.atrous_bdsm_do
            filter = {'tr': {'size': 3, 'vec': [1. / 4, 1. / 2, 1. / 4], 'name': 'Triangle'},
                      'b3': {'size': 5, 'vec': [1. / 16, 1. / 4, 3. / 8, 1. / 4, 1. / 16], 'name': 'B3 spline'}}

            if dobdsm:
                wchain, wopts = self.setpara_bdsm(img)

            n, m = img.ch0_arr.shape

            # Calculate residual image that results from normal (non-wavelet) Gaussian fitting
            Op_make_residimage()(img)
            resid = img.resid_gaus_arr

            lpf = img.opts.atrous_lpf
            if lpf not in ['b3', 'tr']:
                lpf = 'b3'
            jmax = img.opts.atrous_jmax
            l = len(filter[lpf]['vec'])             # 1st 3 is arbit and 2nd 3 is whats expected for a-trous
            if jmax < 1 or jmax > 15:                   # determine jmax
                # Check if largest island size is
                # smaller than 1/3 of image size. If so, use it to determine jmax.
                min_size = min(resid.shape)
                max_isl_shape = (0, 0)
                for isl in img.islands:
                    if isl.image.shape[0] * isl.image.shape[1] > max_isl_shape[0] * max_isl_shape[1]:
                        max_isl_shape = isl.image.shape
                if max_isl_shape != (0, 0) and min(max_isl_shape) < min(resid.shape) / 3.0:
                    min_size = min(max_isl_shape) * 4.0
                else:
                    min_size = min(resid.shape)
                jmax = int(floor(log((min_size / 3.0 * 3.0 - l) / (l - 1) + 1) / log(2.0) + 1.0)) + 1
                if min_size * 0.55 <= (l + (l - 1) * (2 ** (jmax) - 1)):
                    jmax = jmax - 1
            img.wavelet_lpf = lpf
            img.wavelet_jmax = jmax
            mylog.info("Using " + filter[lpf]['name'] + ' filter with J_max = ' + str(jmax))

            img.atrous_islands = []
            img.atrous_gaussians = []
            img.atrous_sources = []
            img.atrous_opts = []
            img.resid_wavelets_arr = cp(img.resid_gaus_arr)

            im_old = img.resid_wavelets_arr
            total_flux = 0.0
            ntot_wvgaus = 0
            stop_wav = False
            pix_masked = N.where(N.isnan(resid))
            jmin = 1
            if img.opts.ncores is None:
                numcores = 1
            else:
                numcores = img.opts.ncores
            for j in range(jmin, jmax + 1):  # extra +1 is so we can do bdsm on cJ as well
                mylogger.userinfo(mylog, "\nWavelet scale #" + str(j))
                im_new = self.atrous(im_old, filter[lpf]['vec'], lpf, j, numcores=numcores, use_scipy_fft=img.opts.use_scipy_fft)
                im_new[pix_masked] = N.nan  # since fftconvolve wont work with blanked pixels
                if img.opts.atrous_sum:
                    w = im_new
                else:
                    w = im_old - im_new
                im_old = im_new
                suffix = 'w' + repr(j)
                filename = img.imagename + '.atrous.' + suffix + '.fits'
                if img.opts.output_all:
                    func.write_image_to_file('fits', filename, w, img, bdir)
                    mylog.info('%s %s' % ('Wrote ', img.imagename + '.atrous.' + suffix + '.fits'))

                # now do bdsm on each wavelet image.
                if dobdsm:
                    wopts['filename'] = filename
                    wopts['basedir'] = bdir
                    box = img.rms_box[0]
                    y1 = (l + (l - 1) * (2 ** (j - 1) - 1))
                    bs = max(5 * y1, box)  # changed from 10 to 5
                    if bs > min(n, m) / 2:
                        wopts['rms_map'] = False
                        wopts['mean_map'] = 'const'
                        wopts['rms_box'] = None
                    else:
                        wopts['rms_box'] = (bs, bs/3)
                        if hasattr(img, '_adapt_rms_isl_pos'):
                            bs_bright = max(5 * y1, img.rms_box_bright[0])
                            if bs_bright < bs/1.5:
                                wopts['adaptive_rms_box'] = True
                                wopts['rms_box_bright'] = (bs_bright, bs_bright/3)
                            else:
                                wopts['adaptive_rms_box'] = False
                    if j <= 3:
                        wopts['ini_gausfit'] = 'default'
                    else:
                        wopts['ini_gausfit'] = 'nobeam'
                    wid = (l + (l - 1) * (2 ** (j - 1) - 1))
                    b1, b2 = img.pixel_beam()[0:2]
                    b1 = b1 * fwsig
                    b2 = b2 * fwsig
                    cdelt = img.wcs_obj.acdelt[:2]

                    wimg = Image(wopts)
                    wimg.beam = (sqrt(wid * wid + b1 * b1) * cdelt[0] * 2.0, sqrt(wid * wid + b2 * b2) * cdelt[1] * 2.0, 0.0)
                    wimg.orig_beam = img.beam
                    wimg.pixel_beam = img.pixel_beam
                    wimg.pixel_beamarea = img.pixel_beamarea
                    wimg.log = 'Wavelet.'
                    wimg.basedir = img.basedir
                    wimg.extraparams['bbsprefix'] = suffix
                    wimg.extraparams['bbsname'] = img.imagename + '.wavelet'
                    wimg.extraparams['bbsappend'] = True
                    wimg.bbspatchnum = img.bbspatchnum
                    wimg.waveletimage = True
                    wimg.j = j
                    wimg.indir = img.indir
                    if hasattr(img, '_adapt_rms_isl_pos'):
                        wimg._adapt_rms_isl_pos = img._adapt_rms_isl_pos

                    self.init_image_simple(wimg, img, w, '.atrous.' + suffix)
                    for op in wchain:
                        op(wimg)
                        gc.collect()
                        if isinstance(op, Op_islands) and img.opts.atrous_orig_isl:
                            if wimg.nisl > 0:

                                # Find islands that do not share any pixels with
                                # islands in original ch0 image.
                                good_isl = []

                                # Make original rank image boolean; rank counts from 0, with -1 being
                                # outside any island
                                orig_rankim_bool = N.array(img.pyrank + 1, dtype=bool)

                                # Multiply rank images
                                old_islands = orig_rankim_bool * (wimg.pyrank + 1) - 1

                                # Exclude islands that don't overlap with a ch0 island.
                                valid_ids = set(old_islands.flatten())
                                for idx, wvisl in enumerate(wimg.islands):
                                    if idx in valid_ids:
                                        wvisl.valid = True
                                        good_isl.append(wvisl)
                                    else:
                                        wvisl.valid = False

                                wimg.islands = good_isl
                                wimg.nisl = len(good_isl)
                                mylogger.userinfo(mylog, "Number of islands found", '%i' % wimg.nisl)

                                # Renumber islands:
                                for wvindx, wvisl in enumerate(wimg.islands):
                                    wvisl.island_id = wvindx

                        if isinstance(op, Op_gausfit):
                            # If opts.atrous_orig_isl then exclude Gaussians outside of
                            # the original ch0 islands
                            nwvgaus = 0
                            if img.opts.atrous_orig_isl:
                                gaul = wimg.gaussians
                                tot_flux = 0.0

                                if img.ngaus == 0:
                                    gaus_id = -1
                                else:
                                    gaus_id = img.gaussians[-1].gaus_num
                                for g in gaul:
                                    if not hasattr(g, 'valid'):
                                        g.valid = False
                                    if not g.valid:
                                        try:
                                            isl_id = img.pyrank[int(g.centre_pix[0] + 1), int(g.centre_pix[1] + 1)]
                                        except IndexError:
                                            isl_id = -1
                                        if isl_id >= 0:
                                            isl = img.islands[isl_id]
                                            gcenter = (int(g.centre_pix[0] - isl.origin[0]),
                                                       int(g.centre_pix[1] - isl.origin[1]))
                                            if not isl.mask_active[gcenter]:
                                                gaus_id += 1
                                                gcp = Gaussian(img, g.parameters[:], isl.island_id, gaus_id)
                                                gcp.gaus_num = gaus_id
                                                gcp.wisland_id = g.island_id
                                                gcp.jlevel = j
                                                g.valid = True
                                                isl.gaul.append(gcp)
                                                isl.ngaus += 1
                                                img.gaussians.append(gcp)
                                                nwvgaus += 1
                                                tot_flux += gcp.total_flux
                                            else:
                                                g.valid = False
                                                g.jlevel = 0
                                        else:
                                            g.valid = False
                                            g.jlevel = 0
                                vg = []
                                for g in wimg.gaussians:
                                    if g.valid:
                                        vg.append(g)
                                wimg.gaussians = vg
                                mylogger.userinfo(mylog, "Number of valid wavelet Gaussians", str(nwvgaus))
                            else:
                                # Keep all Gaussians and merge islands that overlap
                                tot_flux = check_islands_for_overlap(img, wimg)

                                # Now renumber the islands and adjust the rank image before going to next wavelet image
                                renumber_islands(img)

                    total_flux += tot_flux
                    if img.opts.interactive and has_pl:
                        dc = '\033[34;1m'
                        nc = '\033[0m'
                        print(dc + '--> Displaying islands and rms image...' + nc)
                        if max(wimg.ch0_arr.shape) > 4096:
                            print(dc + '--> Image is large. Showing islands only.' + nc)
                            wimg.show_fit(rms_image=False, mean_image=False, ch0_image=False,
                                          ch0_islands=True, gresid_image=False, sresid_image=False,
                                          gmodel_image=False, smodel_image=False, pyramid_srcs=False)
                        else:
                            wimg.show_fit()
                        prompt = dc + "Press enter to continue or 'q' stop fitting wavelet images : " + nc
                        answ = raw_input_no_history(prompt)
                        while answ != '':
                            if answ == 'q':
                                img.wavelet_jmax = j
                                stop_wav = True
                                break
                            answ = raw_input_no_history(prompt)
                    if len(wimg.gaussians) > 0:
                        img.resid_wavelets_arr = self.subtract_wvgaus(img.opts, img.resid_wavelets_arr, wimg.gaussians, wimg.islands)
                        if img.opts.atrous_sum:
                            im_old = self.subtract_wvgaus(img.opts, im_old, wimg.gaussians, wimg.islands)
                    if stop_wav:
                        break

            pyrank = N.zeros(img.pyrank.shape, dtype=N.int32)
            for i, isl in enumerate(img.islands):
                isl.island_id = i
                for g in isl.gaul:
                    g.island_id = i
                for dg in isl.dgaul:
                    dg.island_id = i
                pyrank[tuple(isl.bbox)] += N.invert(isl.mask_active) * (i + 1)
            pyrank -= 1  # align pyrank values with island ids and set regions outside of islands to -1
            img.pyrank = pyrank

            img.ngaus += ntot_wvgaus
            img.total_flux_gaus += total_flux
            mylogger.userinfo(mylog, "Total flux density in model on all scales", '%.3f Jy' % img.total_flux_gaus)
            if img.opts.output_all:
                func.write_image_to_file('fits', img.imagename + '.atrous.cJ.fits',
                                         im_new, img, bdir)
                mylog.info('%s %s' % ('Wrote ', img.imagename + '.atrous.cJ.fits'))
                func.write_image_to_file('fits', img.imagename + '.resid_wavelets.fits',
                                         (img.ch0_arr - img.resid_gaus_arr + img.resid_wavelets_arr), img, bdir + '/residual/')
                mylog.info('%s %s' % ('Wrote ', img.imagename + '.resid_wavelets.fits'))
                func.write_image_to_file('fits', img.imagename + '.model_wavelets.fits',
                                         (img.resid_gaus_arr - img.resid_wavelets_arr), img, bdir + '/model/')
                mylog.info('%s %s' % ('Wrote ', img.imagename + '.model_wavelets.fits'))
            img.completed_Ops.append('wavelet_atrous')

    def atrous(self, image, filtvec, lpf, j, numcores=1, use_scipy_fft=True):

        ff = filtvec[:]
        for i in range(1, len(filtvec)):
            ii = 1 + (2 ** (j - 1)) * (i - 1)
            ff[ii:ii] = [0] * (2 ** (j - 1) - 1)
        kern = N.outer(ff, ff)
        unmasked = N.nan_to_num(image)
        if use_scipy_fft:
            im_new = scipy.signal.fftconvolve(unmasked, kern, mode='same')
        else:
            im_new = fftconvolve(unmasked, kern, mode='same', pad_to_power_of_two=False, numcores=numcores)
        if im_new.shape != image.shape:
            im_new = im_new[0:image.shape[0], 0:image.shape[1]]

        return im_new

    def setpara_bdsm(self, img):
        chain = [Op_preprocess, Op_rmsimage(), Op_threshold(), Op_islands(),
                 Op_gausfit(), Op_gaul2srl(), Op_make_residimage()]

        opts = {'thresh': 'hard'}
        opts['thresh_pix'] = img.thresh_pix
        opts['kappa_clip'] = 3.0
        opts['rms_map'] = img.opts.rms_map
        opts['mean_map'] = img.opts.mean_map
        opts['thresh_isl'] = img.opts.thresh_isl
        opts['minpix_isl'] = 6
        opts['savefits_rmsim'] = False
        opts['savefits_meanim'] = False
        opts['savefits_rankim'] = False
        opts['savefits_normim'] = False
        opts['polarisation_do'] = False
        opts['aperture'] = None
        opts['group_by_isl'] = img.opts.group_by_isl
        opts['quiet'] = img.opts.quiet
        opts['ncores'] = img.opts.ncores

        opts['flag_smallsrc'] = False
        opts['flag_minsnr'] = 0.2
        opts['flag_maxsnr'] = 1.2
        opts['flag_maxsize_isl'] = 2.5
        opts['flag_bordersize'] = 0
        opts['flag_maxsize_bm'] = 50.0
        opts['flag_minsize_bm'] = 0.2
        opts['flag_maxsize_fwhm'] = 0.5
        opts['bbs_patches'] = img.opts.bbs_patches
        opts['filename'] = ''
        opts['output_all'] = img.opts.output_all
        opts['verbose_fitting'] = img.opts.verbose_fitting
        opts['split_isl'] = False
        opts['peak_fit'] = True
        opts['peak_maxsize'] = 30.0
        opts['detection_image'] = ''
        opts['verbose_fitting'] = img.opts.verbose_fitting

        ops = []
        for op in chain:
            if isinstance(op, type):
                ops.append(op())
            else:
                ops.append(op)

        return ops, opts

    def init_image_simple(self, wimg, img, w, name):
        wimg.ch0_arr = w
        wimg.ch0_Q_arr = None
        wimg.ch0_U_arr = None
        wimg.ch0_V_arr = None
        wimg.wcs_obj = img.wcs_obj
        wimg.parentname = img.filename
        wimg.filename = img.filename + name
        wimg.imagename = img.imagename + name + '.pybdsf'
        wimg.pix2sky = img.pix2sky
        wimg.sky2pix = img.sky2pix
        wimg.pix2beam = img.pix2beam
        wimg.beam2pix = img.beam2pix
        wimg.pix2gaus = img.pix2gaus
        wimg.gaus2pix = img.gaus2pix
        wimg.pix2coord = img.pix2coord
        wimg.masked = img.masked
        wimg.mask_arr = img.mask_arr
        wimg.use_io = img.use_io
        wimg.do_cache = img.do_cache
        wimg.tempdir = img.tempdir
        wimg.shape = img.shape
        wimg.frequency = img.frequency
        wimg.equinox = img.equinox
        wimg.use_io = 'fits'

    def subtract_wvgaus(self, opts, residim, gaussians, islands):
        from . import functions as func
        from .make_residimage import Op_make_residimage as opp

        dummy = opp()
        shape = residim.shape
        thresh = opts.fittedimage_clip

        for g in gaussians:
            if g.valid:
                C1, C2 = g.centre_pix
                if hasattr(g, 'wisland_id'):
                    isl = islands[g.wisland_id]
                else:
                    isl = islands[g.island_id]
                b = opp.find_bbox(dummy, thresh * isl.rms, g)
                bbox = N.s_[max(0, int(C1 - b)):min(shape[0], int(C1 + b + 1)),
                            max(0, int(C2 - b)):min(shape[1], int(C2 + b + 1))]
                x_ax, y_ax = N.mgrid[bbox]
                ffimg = func.gaussian_fcn(g, x_ax, y_ax)
                residim[bbox] = residim[bbox] - ffimg

        return residim

    def morphfilter_pyramid(self, img, bdir):
        from math import ceil, floor

        jmax = img.wavelet_jmax
        ind = [i for i, isl in enumerate(img.atrous_islands) if len(isl) > 0]
        ind.reverse()
        lpyr = []
        img.npyrsrc = -1
        if len(ind) > 0:
            for i in ind:
                isls = img.atrous_islands[i]
                for isl in isls:
                    if i != ind[0]:
                        dumr = []
                        for pyrsrc in lpyr:
                            belongs = pyrsrc.belongs(img, isl)
                            if belongs:
                                dumr.append(pyrsrc.pyr_id)
                        if len(dumr) == 1:
                            dumr = dumr[0]
                            pyrsrc = lpyr[dumr]
                            pyrsrc.add_level(img, i, isl)
                        else:
                            pyrsrc = Pyramid_source(img, isl, i)
                            lpyr.append(pyrsrc)
                    else:
                        pyrsrc = Pyramid_source(img, isl, i)
                        lpyr.append(pyrsrc)
        img.pyrsrcs = lpyr

        if img.opts.plot_pyramid and has_pl:
            pl.figure()
            a = ceil(sqrt(jmax))
            b = floor(jmax / a)
            if a * b < jmax:
                b += 1
            colours = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
            sh = img.ch0_arr.shape
            for pyr in img.pyrsrcs:
                for iisl, isl in enumerate(pyr.islands):
                    jj = pyr.jlevels[iisl]
                    col = colours[pyr.pyr_id % 7]
                    pl.subplot(a, b, jj)
                    ind = N.where(~isl.mask_active)
                    pl.plot(ind[0] + isl.origin[0], ind[1] + isl.origin[1], '.', color=col)
                    pl.axis([0.0, sh[0], 0.0, sh[1]])
                    pl.title('J = ' + str(jj))
            pl.savefig(bdir + img.imagename + '.pybdsf.atrous.pyramidsrc.png')


class Pyramid_source(object):
    """ Pyramid_source is a source constructed out of multiple wavelet transform images. """

    def __init__(self, img, island, level0):
        img.npyrsrc = img.npyrsrc + 1
        self.pyr_id = img.npyrsrc
        self.islands = [island]
        self.jlevels = [level0]

    def belongs(self, img, isl):
        from . import functions as func

        # Get centroid of island (as integer)
        mom = func.momanalmask_gaus(isl.image, isl.mask_active, 0, 1.0, False)
        cen = N.array(mom[1:3]) + isl.origin
        belong = False

        # Check if lies within any island of self
        for i, pyrisl in enumerate(self.islands):
            if N.sum([pyrisl.bbox[j].start <= cen[j] < pyrisl.bbox[j].stop for j in range(2)]) == 2:
                pix = tuple([cen[j] - pyrisl.origin[j] for j in range(2)])
                if not pyrisl.mask_active[pix]:
                    belong = True

        return belong

    def add_level(self, img, level, isl):
        self.islands.append(isl)
        self.jlevels.append(level + 1)


Image.pyrsrcs = List(tInstance(Pyramid_source), doc="List of Pyramidal sources")


def fftconvolve(in1, in2, mode="full", pad_to_power_of_two=True, numcores=1):
    """Convolve two N-dimensional arrays using FFT. See convolve.

    """
    s1 = array(in1.shape)
    s2 = array(in2.shape)
    complex_result = (N.issubdtype(in1.dtype, N.complex) or
                      N.issubdtype(in2.dtype, N.complex))
    size = s1 + s2 - 1

    if pad_to_power_of_two:
        # Use 2**n-sized FFT; it might improve performance
        fsize = 2 ** N.ceil(N.log2(size))
    else:
        # Padding to a power of two might degrade performance, too
        fsize = size
    if has_pyfftw:
        IN1 = N.fft.fftn(in1, fsize, threads=numcores)
        IN1 *= N.fft.fftn(in2, fsize, threads=numcores)
        fslice = tuple([slice(0, int(sz)) for sz in size])
        ret = N.fft.ifftn(IN1, threads=numcores)[fslice].copy()
    else:
        IN1 = N.fft.fftn(in1, fsize)
        IN1 *= N.fft.fftn(in2, fsize)
        fslice = tuple([slice(0, int(sz)) for sz in size])
        ret = N.fft.ifftn(IN1)[fslice].copy()
    del IN1
    if not complex_result:
        ret = ret.real
    if mode == "full":
        return ret
    elif mode == "same":
        if product(s1, axis=0) > product(s2, axis=0):
            osize = s1
        else:
            osize = s2
        return func.centered(ret, osize)
    elif mode == "valid":
        return func.centered(ret, abs(s2 - s1) + 1)


def rebase_bbox(box, minv):
    # return a new bounding box tuple where minv is subtracted from
    # all the co-ordinate values
    nbox = []
    for i, sl in enumerate(box):
        nbox.append(slice(sl.start-minv[i], sl.stop-minv[i], None))
    return tuple(nbox)


def merge_bbox(box1, box2):
    # For two bounding box tuples find the minimal n-dimensional space
    # that encompasses both structures and make new bounding boxes in
    # this co-ordinate system
    minv = []
    maxv = []
    for sl1, sl2 in zip(box1, box2):
        minv.append(min(sl1.start, sl2.start))
        maxv.append(max(sl1.stop, sl2.stop))
    nbox1 = rebase_bbox(box1, minv)
    nbox2 = rebase_bbox(box2, minv)
    dims = [y-x for x, y in zip(minv, maxv)]
    fullbox = [slice(x, y, None) for x, y in zip(minv, maxv)]
    return dims, nbox1, nbox2, N.array(minv), fullbox


def merge_islands(img, isl1, isl2):
    """Merge two islands into one

    Final island has island_id of isl1. The Gaussians from isl2 are appended
    those in the isl1 list, with numbering starting from the last number in
    img.gaussians (which is also updated with the isl2 Gaussians).

    The merged island replaces isl1 in img.
    """
    from .islands import Island
    import scipy.ndimage as nd

    shape, nbox1, nbox2, origin, fullbox = merge_bbox(isl1.bbox, isl2.bbox)
    mask1 = N.zeros(shape, dtype=bool)
    mask1[nbox1] = ~isl1.mask_active
    mask2 = N.zeros(shape, dtype=bool)
    mask2[nbox2] = ~isl2.mask_active
    overlap_mask = N.logical_and(mask1, mask2)
    if N.any(overlap_mask):
        full_mask = N.logical_or(mask1, mask2)
        image = img.ch0_arr
        mask = img.mask_arr
        rms = img.rms_arr
        mean = img.mean_arr
        rank = len(image.shape)
        connectivity = nd.generate_binary_structure(rank, rank)
        labels, count = nd.label(full_mask, connectivity)
        slices = nd.find_objects(labels)
        bbox = slices[0]
        new_bbox = rebase_bbox(bbox, -origin)
        idx = isl1.island_id
        # labels array passed to Island must be capable of being
        # indexed by new bounding box, so convert. Do the subtraction
        # first to avoid an expensive operation over the whole array
        labels = labels-1+idx
        new_labels = N.zeros(image.shape)
        new_labels[tuple(fullbox)] = labels
        beamarea = img.pixel_beamarea()
        merged_isl = Island(image, mask, mean, rms, new_labels, new_bbox, idx, beamarea)

        # Add all the Gaussians to the merged island
        merged_isl.gaul = isl1.gaul
        merged_isl.dgaul = isl1.dgaul
        copy_gaussians(img, merged_isl, isl2)
        img.islands[idx] = merged_isl


def copy_gaussians(img, isl1, isl2):
    """Copies Gaussians from isl2 to isl1

    img.gaussians is also updated
    """
    if img.ngaus == 0:
        gaus_id = -1
    else:
        gaus_id = img.gaussians[-1].gaus_num
    for g in isl2.gaul:
        gaus_id += 1
        gcp = Gaussian(img, g.parameters[:], isl1.island_id, gaus_id)
        gcp.gaus_num = gaus_id
        gcp.jlevel = g.jlevel
        if g.jlevel > 0:
            # Preserve the wavelet rms and mean values if the isl2 Gaussian was fit to
            # a wavelet image
            gcp.wave_rms = g.rms
            gcp.wave_mean = g.mean
        else:
            gcp.wave_rms = 0.0
            gcp.wave_mean = 0.0

        isl1.gaul.append(gcp)
        img.ngaus += 1
        img.gaussians.append(gcp)


def renumber_islands(img):
    """Renumbers island_ids (after, e.g., removing one)

    Also renumbers the pyrank image.
    """
    pyrank = N.zeros(img.pyrank.shape, dtype=N.int32)
    for i, isl in enumerate(img.islands):
        isl.island_id = i
        for g in isl.gaul:
            g.island_id = i
        for dg in isl.dgaul:
            dg.island_id = i
            pyrank[tuple(isl.bbox)] += N.invert(isl.mask_active) * (i + 1)
    pyrank -= 1  # align pyrank values with island ids and set regions outside of islands to -1
    img.pyrank = pyrank
    gaussian_list = [g for isl in img.islands for g in isl.gaul]
    img.gaussians = gaussian_list


def check_islands_for_overlap(img, wimg):

    """Checks for overlaps between img and wimg islands"""

    have_numexpr = True
    try:
        import numexpr as ne
    except:
        have_numexpr = False

    tot_flux = 0.0
    bar = statusbar.StatusBar('Checking islands for overlap ............ : ', 0, len(wimg.islands))

    # Make masks for regions that have islands
    wpp = wimg.pyrank+1  # does not change, store for later
    wav_rankim_bool = wpp > 0  # boolean
    orig_rankim_bool = img.pyrank > -1

    # Make "images" of island ids for overlaping regions
    orig_islands = wav_rankim_bool * (img.pyrank + 1) - 1
    if not img.opts.quiet:
        bar.start()
    for idx, wvisl in enumerate(wimg.islands):
        if len(wvisl.gaul) > 0:
            # Get unique island IDs. If an island overlaps with one
            # in the original ch0 image, merge them together. If not,
            # add the island as a new one.
            wav_islands = orig_rankim_bool[tuple(wvisl.bbox)] * wpp[tuple(wvisl.bbox)] - 1
            wav_ids = N.unique(wav_islands)  # saves conversion to set and back
            for wvg in wvisl.gaul:
                tot_flux += wvg.total_flux
                wvg.valid = True
            if idx in wav_ids:
                orig_idx = N.unique(orig_islands[tuple(wvisl.bbox)][wav_islands == idx])
                if len(orig_idx) == 1:
                    merge_islands(img, img.islands[orig_idx[0]], wvisl)
                else:
                    merge_islands(img, img.islands[orig_idx[0]], wvisl)
                    for oidx in orig_idx[1:]:
                        merge_islands(img, img.islands[orig_idx[0]], img.islands[oidx])
                    img.islands = [x for x in img.islands if x.island_id not in orig_idx[1:]]
                    renumber_islands(img)
                # Now recalculate the overlap images, since the islands have changed
                ipp = img.pyrank+1
                if have_numexpr:
                    orig_islands = ne.evaluate('wav_rankim_bool * ipp - 1')
                else:
                    orig_islands = wav_rankim_bool * ipp - 1
            else:
                isl_id = img.islands[-1].island_id + 1
                new_isl = wvisl.copy(img.pixel_beamarea(), image=img.ch0_arr[tuple(wvisl.bbox)],
                                     mean=img.mean_arr[tuple(wvisl.bbox)],
                                     rms=img.rms_arr[tuple(wvisl.bbox)])
                new_isl.gaul = []
                new_isl.dgaul = []
                new_isl.island_id = isl_id
                img.islands.append(new_isl)
                copy_gaussians(img, new_isl, wvisl)

        if not img.opts.quiet:
            bar.increment()
    bar.stop()

    return tot_flux
