"""Module islands.

Defines operation Op_islands which does island detection.
Current implementation uses scipy.ndimage operations for island detection.
While it's implemented to work for images of arbitrary dimensionality,
the bug in the current version of scipy (0.6) often causes crashes
(or just wrong results) for 3D inputs.

If this (scipy.ndimage.label) isn't fixed by the time we need 3D source
extraction, one will have to adopt my old pixel-runs algorithm for 3D data.
Check out islands.py rev. 1362 from repository for it.
"""
from __future__ import absolute_import
from __future__ import division

import numpy as N
import scipy.ndimage as nd
from .image import *
from . import mylogger
from . import functions as func
from .output import write_islands
from .readimage import Op_readimage
from .preprocess import Op_preprocess
from .rmsimage import Op_rmsimage
from .threshold import Op_threshold
from .collapse import Op_collapse
import os


class Op_islands(Op):
    """Detect islands of emission in the image

    All detected islands are stored in the list img.islands,
    where each individual island is represented as an instance
    of class Island.

    The option to detect islands on a different "detection"
    image is also available. This option is useful for example
    when a primary beam correction is used -- it is generally
    better to detect sources on the uncorrected image, but
    to measure them on the corrected image.

    Prerequisites: module rmsimage should be run first.
    """
    def __call__(self, img):
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Islands")
        opts = img.opts

        minsize = opts.minpix_isl
        if minsize is None:
            minsize = int(img.pixel_beamarea()/3.0)  # 1/3 of beam area in pixels
            if minsize < 6:
                minsize = 6  # Need at least 6 pixels to obtain good fits
            mylogger.userinfo(mylog, "Minimum number of pixels per island", '%i' %
                              minsize)
        img.minpix_isl = minsize
        maxsize = opts.maxpix_isl
        if maxsize is None:
            maxsize = N.inf
        img.maxpix_isl = maxsize

        if opts.detection_image != '':
            # Use a different image for island detection. The detection
            # image and the measurement image must have the same shape
            # and be registered. Otherwise, one could reproject the
            # detection image using, e.g., the Kapteyn package.
            #
            # First, set up up an Image object and run a limited
            # op_chain.
            from . import _run_op_list
            mylogger.userinfo(mylog, "\nDetermining islands from detection image")

            det_chain, det_opts = self.setpara_bdsm(img, opts.detection_image)
            det_img = Image(det_opts)
            det_img.log = 'Detection image'
            success = _run_op_list(det_img, det_chain)
            if not success:
                return

            # Check that the ch0 images are the same size
            ch0_map = img.ch0_arr
            det_ch0_map = det_img.ch0_arr
            det_shape = det_ch0_map.shape
            ch0_shape = ch0_map.shape
            if det_shape != ch0_shape:
                raise RuntimeError("Detection image shape does not match that of input image.")

            # Save the rms and mean maps derived from the detection image
            img.detection_mean_arr = det_img.mean_arr
            img.detection_rms_arr = det_img.rms_arr

            # Run through islands and correct the image and rms, mean and max values
            corr_islands = []
            mean_map = img.mean_arr
            rms_map = img.rms_arr
            for i, isl in enumerate(det_img.islands):
                islcp = isl.copy(img.pixel_beamarea(), image=ch0_map[tuple(isl.bbox)], mean=mean_map[tuple(isl.bbox)], rms=rms_map[tuple(isl.bbox)])
                islcp.island_id = i
                corr_islands.append(islcp)
            img.islands = corr_islands
            img.nisl = len(img.islands)
            img.pyrank = det_img.pyrank
            img.minpix_isl = det_img.minpix_isl

            if opts.savefits_det_rmsim or opts.output_all:
                resdir = os.path.join(img.basedir, 'background')
                os.makedirs(resdir, exist_ok=True)
                func.write_image_to_file(img.use_io, img.imagename + '.detection_rmsd_I.fits',
                                         img.detection_rms_arr, img, resdir)
                mylog.info('%s %s' % ('Writing ', os.path.join(resdir, img.imagename+'.detection_rmsd_I.fits')))
            if opts.savefits_det_meanim or opts.output_all:
                resdir = os.path.join(img.basedir, 'background')
                os.makedirs(resdir, exist_ok=True)
                func.write_image_to_file(img.use_io, img.imagename + '.detection_mean_I.fits',
                                         img.detection_mean_arr, img, resdir)
                mylog.info('%s %s' % ('Writing ', os.path.join(resdir, img.imagename+'.detection_mean_I.fits')))

            mylogger.userinfo(mylog, "\nContinuing processing using primary image")
        else:
            if opts.src_ra_dec is not None:
                mylogger.userinfo(mylog, "Constructing islands at user-supplied source locations")
                img.islands = self.coords_to_isl(img, opts)
            else:
                img.islands = self.ndimage_alg(img, opts)
            img.nisl = len(img.islands)

            mylogger.userinfo(mylog, "Number of islands found", '%i' %
                              len(img.islands))

            ch0_map = img.ch0_arr
            ch0_shape = ch0_map.shape
            pyrank = N.zeros(ch0_shape, dtype=N.int32)
            for i, isl in enumerate(img.islands):
                isl.island_id = i
                pyrank[tuple(isl.bbox)] += N.invert(isl.mask_active) * (i + 1)
            pyrank -= 1  # align pyrank values with island ids and set regions outside of islands to -1

            img.pyrank = pyrank

        if opts.output_all:
            write_islands(img)
        if opts.savefits_rankim or opts.output_all:
            resdir = os.path.join(img.basedir, 'misc')
            os.makedirs(resdir, exist_ok=True)
            func.write_image_to_file(img.use_io, img.imagename + '_pyrank.fits', img.pyrank, img, resdir)
            mylog.info('%s %s' % ('Writing ', os.path.join(resdir, img.imagename+'_pyrank.fits')))

        img.completed_Ops.append('islands')
        return img

    def ndimage_alg(self, img, opts):
        """Island detection using scipy.ndimage

        Use scipy.ndimage.label to detect islands of emission in the image.
        Island is defined as group of tightly connected (8-connectivity
        for 2D images) pixels with emission.

        The following cuts are applied:
         - pixel is considered to have emission if it is 'thresh_isl' times
           higher than RMS.
         - Island should have at least 'minsize' active pixels
         - There should be at lease 1 pixel in the island which is 'thresh_pix'
           times higher than noise (peak clip).

        Parameters:
        image, mask: arrays with image data and mask
        mean, rms: arrays with mean & rms maps
        thresh_isl: threshold for 'active pixels'
        thresh_pix: threshold for peak
        minsize: minimal acceptable island size

        Function returns a list of Island objects.
        """
        # Islands detection
        image = img.ch0_arr
        mask = img.mask_arr
        rms = img.rms_arr
        mean = img.mean_arr
        thresh_isl = opts.thresh_isl
        thresh_pix = img.thresh_pix

        # Here act_pixels is true if significant emission
        if img.masked:
            act_pixels = ~(mask.copy())
            act_pixels[~mask] = (image[~mask]-mean[~mask])/thresh_isl >= rms[~mask]
        else:
            act_pixels = (image-mean)/thresh_isl >= rms

        # Find dimension of image
        rank = len(image.shape)
        # Generates matrix for connectivity, in this case, 8-conn
        connectivity = nd.generate_binary_structure(rank, rank)
        # Here labels = matrix with value = (initial) island number
        labels, count = nd.label(act_pixels, connectivity)
        # Here slices has limits of bounding box of each such island
        slices = nd.find_objects(labels)
        img.island_labels = labels

        # Apply cuts on island size and peak value
        pyrank = N.zeros(image.shape, dtype=N.int32)
        res = []
        for idx, s in enumerate(slices):
            idx += 1  # nd.labels indices are counted from 1
            isl_size = (labels[s] == idx).sum()  # number of pixels inside bounding box which are in island
            isl_peak = nd.maximum(image[s], labels[s], idx)
            isl_maxposn = tuple(N.array(N.unravel_index(N.nanargmax(image[s]), image[s].shape)) +
                                N.array((s[0].start, s[1].start)))
            if (isl_size >= img.minpix_isl) and (isl_size <= img.maxpix_isl) and (isl_peak - mean[isl_maxposn])/thresh_pix > rms[isl_maxposn]:
                isl = Island(image, mask, mean, rms, labels, s, idx, img.pixel_beamarea())
                res.append(isl)
                pyrank[tuple(isl.bbox)] += N.invert(isl.mask_active)*idx // idx

        return res

    def coords_to_isl(self, img, opts):
        """Construct islands around given coordinates with given size.

        Returns a list of island objects.
        """
        coords = opts.src_ra_dec  # list of RA and Dec tuples
        isl_radius_pix = opts.src_radius_pix
        if isl_radius_pix is None:
            isl_radius_pix = img.beam2pix(img.beam)[0]  # twice beam major axis radius at half max (= FWHM)

        res = []
        for idx, coord in enumerate(coords):
            idx += 1  # nd.labels indices are counted from 1
            isl_posn_pix = img.sky2pix(coord)
            image = img.ch0_arr
            mask = img.mask_arr
            rms = img.rms_arr
            mean = img.mean_arr
            labels = func.make_src_mask(image.shape, isl_posn_pix, isl_radius_pix)
            if img.masked:
                aper_mask = N.where(labels.astype(bool) & ~mask)
            else:
                aper_mask = N.where(labels.astype(bool))
            if N.size(aper_mask) >= img.minpix_isl and N.size(aper_mask) <= img.maxpix_isl:
                labels[aper_mask] = idx
                s = [slice(max(0, isl_posn_pix[0] - isl_radius_pix - 1),
                     min(image.shape[0], isl_posn_pix[0] + isl_radius_pix + 1)),
                     slice(max(0, isl_posn_pix[1] - isl_radius_pix - 1),
                     min(image.shape[1], isl_posn_pix[1] + isl_radius_pix + 1))]
                isl = Island(image, mask, mean, rms, labels, s, idx, img.pixel_beamarea())
                res.append(isl)
        return res

    def setpara_bdsm(self, img, det_file):

        chain = [Op_readimage(), Op_collapse(), Op_preprocess, Op_rmsimage(),
                 Op_threshold(), Op_islands()]
        opts = img.opts.to_dict()
        opts['filename'] = det_file
        opts['detection_image'] = ''
        opts['polarisation_do'] = False
        opts['rmsmean_map_filename'] = opts['rmsmean_map_filename_det']
        opts['det_rmsmean_map_filename'] = None

        ops = []
        for op in chain:
            if isinstance(op, type):
                ops.append(op())
            else:
                ops.append(op)

        return ops, opts


class Island(object):
    """Instances of this class represent islands of emission in the image.

    Its primary use is a container for all kinds of data describing island.
    """
    def __init__(self, img, mask, mean, rms, labels, bbox, idx,
                 beamarea, origin=None, noise_mask=None, copy=False):
        """Create Island instance.

        Parameters:
        img, mask, mean, rms: arrays describing image
        labels: labels array from scipy.ndimage
        bbox: slices
        """
        # Add attribute definitions needed for output
        self.island_id_def = Int(doc="Island id, starting from 0", colname='Isl_id')
        self.shapelet_basis_def = String(doc="Coordinate system for shapelet decomposition (cartesian/polar)",
                                         colname='Basis', units=None)
        self.shapelet_beta_def = Float(doc="Value of shapelet scale beta", colname='Beta', units=None)
        self.shapelet_nmax_def = Int(doc="Maximum value of shapelet order", colname='NMax', units=None)
        self.shapelet_posn_sky_def = List(Float(), doc="Posn (RA, Dec in deg) of shapelet centre",
                                          colname=['RA', 'DEC'], units=['deg', 'deg'])
        self.shapelet_posn_skyE_def = List(Float(), doc="Error on sky coordinates of shapelet centre",
                                           colname=['E_RA', 'E_DEC'], units=['deg', 'deg'])
        self.shapelet_cf_def = NArray(doc="Coefficient matrix of the shapelet decomposition",
                                      colname='Coeff_matrix', units=None)

        if not copy:
            # We make bbox slightly bigger
            self.oldbbox = bbox
            self.oldidx = idx
            bbox = self.__expand_bbox(bbox, img.shape)
            origin = [b.start for b in bbox]   # easier in case ndim > 2
            data = img[tuple(bbox)]
            bbox_rms_im = rms[tuple(bbox)]
            bbox_mean_im = mean[tuple(bbox)]

            # Create (inverted) masks
            # Note that mask_active is the island mask; mask_noisy marks only
            # the noisy pixels in the island image. If you want to mask the
            # noisy pixels, set the final mask to:
            #     mask = mask_active + mask_noisy
            isl_mask = (labels[tuple(bbox)] == idx)
            noise_mask = (labels[tuple(bbox)] == 0)
            N.logical_or(noise_mask, isl_mask, noise_mask)

            # Invert masks
            N.logical_not(isl_mask, isl_mask)
            N.logical_not(noise_mask, noise_mask)
            if isinstance(mask, N.ndarray):
                noise_mask[mask[tuple(bbox)]] = True
                isl_mask[mask[tuple(bbox)]] = True
        else:
            if origin is None:
                origin = [b.start for b in bbox]
            isl_mask = mask
            if noise_mask is None:
                noise_mask = mask
            data = img
            bbox_rms_im = rms
            bbox_mean_im = mean
            self.oldbbox = bbox
            self.oldidx = idx

        # Finish initialization
        isl_size = N.sum(~isl_mask)
        self.island_id = idx
        self.bbox = bbox
        self.origin = origin
        self.image = data
        self.mask_active = isl_mask
        self.mask_noisy = noise_mask
        self.shape = data.shape
        self.size_active = isl_size
        self.max_value = N.max(self.image[~self.mask_active])
        in_bbox_and_unmasked = N.where(~N.isnan(bbox_rms_im))
        self.rms = bbox_rms_im[in_bbox_and_unmasked].mean()
        in_bbox_and_unmasked = N.where(~N.isnan(bbox_mean_im))
        self.mean = bbox_mean_im[in_bbox_and_unmasked].mean()
        self.islmean = bbox_mean_im[in_bbox_and_unmasked].mean()
        self.total_flux = N.nansum(self.image[in_bbox_and_unmasked])/beamarea
        pixels_in_isl = N.sum(~N.isnan(self.image[self.mask_active]))  # number of unmasked pixels assigned to current island
        self.total_fluxE = func.nanmean(bbox_rms_im[in_bbox_and_unmasked]) * N.sqrt(pixels_in_isl/beamarea)  # Jy
        self.border = self.get_border()
        self.gaul = []
        self.fgaul = []
        self.sources = []
        self.gresid_mean = 0.0
        self.gresid_rms = 0.0

    def __setstate__(self, state):
        """Needed for multiprocessing"""
        self.mean = state['mean']
        self.rms = state['rms']
        self.image = state['image']
        self.islmean = state['islmean']
        self.mask_active = state['mask_active']
        self.mask_noisy = state['mask_noisy']
        self.size_active = state['size_active']
        self.shape = state['shape']
        self.origin = state['origin']
        self.island_id = state['island_id']
        self.oldidx = state['oldidx']
        self.bbox = state['bbox']

    def __getstate__(self):
        """Needed for multiprocessing"""
        state = {}
        state['mean'] = self.mean
        state['rms'] = self.rms
        state['image'] = self.image
        state['islmean'] = self.islmean
        state['mask_active'] = self.mask_active
        state['mask_noisy'] = self.mask_noisy
        state['size_active'] = self.size_active
        state['shape'] = self.shape
        state['origin'] = self.origin
        state['island_id'] = self.island_id
        state['oldidx'] = self.oldidx
        state['bbox'] = self.bbox
        return state

    # Do map etc in case of ndim image
    def __expand_bbox(self, bbox, shape):
        """Expand bbox of the image by 1 pixel"""
        def __expand(bbox, shape):
            return slice(int(max(0, bbox.start - 1)), int(min(shape, bbox.stop + 1)))
        ebbox = [__expand(b, shape[i]) for i, b in enumerate(bbox)]
        return ebbox

    def copy(self, pixel_beamarea, image=None, mean=None, rms=None):
        mask = self.mask_active
        noise_mask = self.mask_noisy
        if image is None:
            image = self.image
        if mean is None:
            mean = N.zeros(mask.shape, dtype=N.float32) + self.mean
        if rms is None:
            rms = N.zeros(mask.shape, dtype=N.float32) + self.rms

        bbox = self.bbox
        idx = self.oldidx
        origin = self.origin
        return Island(image, mask, mean, rms, None, bbox, idx, pixel_beamarea,
                      origin=origin, noise_mask=noise_mask, copy=True)

    def get_border(self):
        """ From all valid island pixels, generate the border."""
        mask = ~self.mask_active
        border = N.transpose(N.asarray(N.where(mask ^ nd.binary_erosion(mask)))) + self.origin

        return N.transpose(N.array(border))
