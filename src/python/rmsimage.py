"""Module rmsimage.

Defines operation Op_rmsimage which calculates mean and
rms maps.

The current implementation will handle both 2D and 3D images,
where for 3D case it will calculate maps for each plane (=
Stokes images).
"""

import numpy as N
import scipy.ndimage as nd
import _cbdsm
from image import Op, Image, NArray, List
import const
import mylogger
import os
import functions as func
import scipy.ndimage as nd
import multi_proc as mp
import itertools

### insert into Image tc-variables for mean & rms maps
Image.mean = NArray(doc="Mean map, Stokes I")
Image.rms  = NArray(doc="RMS map, Stokes I")
Image.mean_QUV = List(NArray(), doc="Mean maps, Stokes QUV")
Image.rms_QUV  = List(NArray(), doc="RMS maps, Stokes QUV")

class Op_rmsimage(Op):
    """Calculate rms & noise maps

    Prerequisites: Module preprocess should be run first.
    """
    def __call__(self, img):
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"RMSimage")
        mylogger.userinfo(mylog, "Calculating background rms and mean images")
        if img.opts.polarisation_do:
            pols = ['I', 'Q', 'U', 'V']
        else:
            pols = ['I'] # assume I is always present

        if hasattr(img, 'rms_mask'):
            mask = img.rms_mask
        else:
            mask = img.mask
        opts = img.opts
        ch0_images = [img.ch0, img.ch0_Q, img.ch0_U, img.ch0_V]
        cmeans = [img.clipped_mean] + img.clipped_mean_QUV
        crmss = [img.clipped_rms] + img.clipped_rms_QUV
        cdelt = N.array(img.wcs_obj.acdelt[:2])

        # Determine box size for rms/mean map calculations.
        # If user specifies rms_box, use it. Otherwise, use either an
        # adaptive binning scheme that shrinks the box near
        # the brightest sources or estimate rms_box from bright sources.
        #
        # The adaptive scheme calculates the rms/mean map
        # at two different scales:
        #   1) using a large rms_box, set by size of largest source
        #   2) using a small rms_box, set by size of largest bright source
        # Then, the rms and mean values at a given point are determined
        # by a weighted average of the values in the maps at the two
        # scales.
        fwsig = const.fwsig
        min_adapt_threshold = 10.0
        if opts.adaptive_thresh == None:
            adapt_thresh = 50.0
            start_thresh = 500.0
        else:
            adapt_thresh = opts.adaptive_thresh
            if adapt_thresh < min_adapt_threshold:
                adapt_thresh = min_adapt_threshold
                opts.adaptive_thresh = min_adapt_threshold
            start_thresh = adapt_thresh
        brightsize = None
        isl_pos = []
        do_adapt = img.opts.adaptive_rms_box
        img.use_rms_map = None
        img.mean_map_type = None

        # 'size' of brightest source
        kappa1 = 3.0
        try:
            brightsize = int(round(2.*img.beam[0]/cdelt[0]/fwsig*
                               sqrt(2.*log(img.max_value/(kappa1*crms)))))
        except:
            brightsize = int(round(2.*img.beam[0]/cdelt[0]/fwsig))
        mylog.info('Estimated size of brightest source (pixels) = '+str(brightsize))

        # Using clipped mean and rms and a starting threshold of 500 sigma,
        # search for bright sources. If fewer than 5 are found, reduce
        # threshold until limit set by adapt_thresh is hit.
        cmean = cmeans[0]
        crms = crmss[0]
        image = ch0_images[0]
        shape = image.shape
        isl_size_bright = []
        isl_area_highthresh = []
        isl_peak = []
        max_isl_brightsize = 0.0
        threshold = start_thresh
        if do_adapt:
            mylogger.userinfo(mylog, "Using adaptive scaling of rms_box")
            while len(isl_size_bright) < 5 and threshold >= adapt_thresh:
                isl_size_bright=[]
                isl_maxposn = []
                act_pixels = (img.ch0-cmean)/threshold >= crms
                threshold *= 0.8
                if isinstance(mask, N.ndarray):
                    act_pixels[mask] = False
                rank = len(image.shape)
                connectivity = nd.generate_binary_structure(rank, rank)
                labels, count = nd.label(act_pixels, connectivity)
                slices = nd.find_objects(labels)
                for idx, s in enumerate(slices):
                    isl_size_bright.append(max([s[0].stop-s[0].start, s[1].stop-s[1].start]))
                    size_area = (labels[s] == idx+1).sum()/img.pixel_beamarea*2.0
                    isl_area_highthresh.append(size_area)
                    isl_maxposn.append(tuple(N.array(N.unravel_index(N.argmax(image[s]), image[s].shape))+\
                          N.array((s[0].start, s[1].start))))
                    isl_peak.append(nd.maximum(image[s], labels[s], idx+1))

        # Check islands found above at thresh_isl threshold to determine if
        # the bright source is embedded inside a large island or not. If it is,
        # exclude it from the bright-island list. Also find the size of the
        # largest island at this threshold to set the large-scale rms_box
        bright_threshold = threshold
        threshold = 10.0
        act_pixels = (img.ch0-cmean)/threshold >= crms
        if isinstance(mask, N.ndarray):
            act_pixels[mask] = False
        rank = len(image.shape)
        connectivity = nd.generate_binary_structure(rank, rank)
        labels, count = nd.label(act_pixels, connectivity)
        slices = nd.find_objects(labels)
        isl_size = []
        isl_size_highthresh = []
        isl_size_lowthresh = []
        isl_snr = []
        thratio = threshold/bright_threshold
        for idx, s in enumerate(slices):
            isl_area_lowthresh = (labels[s] == idx+1).sum()/img.pixel_beamarea*2.0
            isl_maxposn_lowthresh = tuple(N.array(N.unravel_index(N.argmax(image[s]), image[s].shape))+
                                          N.array((s[0].start, s[1].start)))
            isl_size += [s[0].stop-s[0].start, s[1].stop-s[1].start]
            if do_adapt and isl_maxposn_lowthresh in isl_maxposn:
                bright_indx = isl_maxposn.index(isl_maxposn_lowthresh)
                if isl_area_lowthresh < 25.0 or isl_area_lowthresh/isl_area_highthresh[bright_indx] < 8.0:
                    isl_pos.append(isl_maxposn_lowthresh)
                    isl_size_lowthresh.append(max([s[0].stop-s[0].start, s[1].stop-s[1].start]))
                    isl_size_highthresh.append(isl_size_bright[bright_indx])
                    isl_snr.append(isl_peak[bright_indx]/crms)

        if len(isl_size) == 0:
            max_isl_size = 0.0
        else:
            max_isl_size = max(isl_size)
        mylog.info('Maximum extent of largest 10-sigma island using clipped rms (pixels) = '+str(max_isl_size))
        if len(isl_size_highthresh) == 0:
            max_isl_size_highthresh = 0.0
            max_isl_size_lowthresh = 0.0
        else:
            max_isl_size_highthresh = max(isl_size_highthresh)
            max_isl_size_lowthresh = max(isl_size_lowthresh)
            avg_max_isl_size = (max_isl_size_highthresh + max_isl_size_lowthresh) / 2.0

        if len(isl_pos) == 0:
            # No bright sources found
            do_adapt = False
        min_size_allowed = int(img.pixel_beam[0]*9.0)

        if opts.rms_box is None or (opts.rms_box_bright is None and do_adapt):
            if do_adapt:
                bsize = int(max(brightsize, min_size_allowed, max_isl_size_highthresh*2.0))
            else:
                bsize = int(max(brightsize, min_size_allowed, max_isl_size*2.0))
            bsize2 = int(max(min(img.ch0.shape)/10.0, max_isl_size*5.0))
            if bsize < min_size_allowed:
                bsize = min_size_allowed
            if bsize % 10 == 0: bsize += 1
            if bsize2 < min_size_allowed:
                bsize2 = min_size_allowed
            if bsize2 % 10 == 0: bsize2 += 1
            bstep = int(round(min(bsize/3., min(shape)/10.)))
            bstep2 = int(round(min(bsize2/3., min(shape)/10.)))
            if opts.rms_box_bright is None:
                img.rms_box_bright = (bsize, bstep)
            else:
                img.rms_box_bright = opts.rms_box_bright
            if opts.rms_box is None:
                img.rms_box = (bsize2, bstep2)
            else:
                img.rms_box = opts.rms_box
        else:
            if do_adapt:
                img.rms_box_bright = opts.rms_box_bright
                img.rms_box = opts.rms_box
            else:
                img.rms_box_bright = opts.rms_box
                img.rms_box = opts.rms_box
        if do_adapt:
            map_opts = (opts.kappa_clip, img.rms_box_bright, opts.spline_rank)
        else:
            map_opts = (opts.kappa_clip, img.rms_box, opts.spline_rank)

        for ipol, pol in enumerate(pols):
          data = ch0_images[ipol]
          mean = N.zeros(data.shape, dtype=N.float32)
          rms  = N.zeros(data.shape, dtype=N.float32)
          if len(pols) > 1:
              pol_txt = ' (' + pol + ')'
          else:
              pol_txt = ''

          ## calculate rms/mean maps if needed
          if ((opts.rms_map is not False) or (opts.mean_map not in ['zero', 'const'])) and img.rms_box[0] > min(img.ch0.shape)/2.0:
            # rms box is too large - just use constant rms and mean
            self.output_rmsbox_size(img)
            mylogger.userinfo(mylog, 'Size of rms_box larger than 1/2 of image size')
            mylogger.userinfo(mylog, 'Using constant background rms and mean')
            img.use_rms_map = False
            img.mean_map_type = 'const'
          else:
            if (opts.rms_map is not False) or (opts.mean_map not in ['zero', 'const']):
              if len(data.shape) == 2:   ## 2d case
                mean, rms = self.calculate_maps(img, data, mean, rms, mask, map_opts, do_adapt=do_adapt,
                                bright_pt_coords=isl_pos, rms_box2=img.rms_box,
                                logname="PyBDSM."+img.log, ncores=img.opts.ncores)
              elif len(data.shape) == 3: ## 3d case
                if not isinstance(mask, N.ndarray):
                  mask = N.zeros(data.shape[0], dtype=bool)
                for i in range(data.shape[0]):
                    ## iterate each plane
                    mean, rms = self.calculate_maps(img, data[i], mean[i], rms[i], mask[i], map_opts,
                                    do_adapt=do_adapt, bright_pt_coords=isl_pos,
                                    rms_box2=img.rms_box, logname="PyBDSM."+img.log,
                                    ncores=img.opts.ncores)
              else:
                mylog.critical('Image shape not handleable' + pol_txt)
                raise RuntimeError("Can't handle array of this shape" + pol_txt)
              self.output_rmsbox_size(img)
              if do_adapt:
                  mylogger.userinfo(mylog, 'Number of sources using small scale', str(len(isl_pos)))
              mylog.info('Background rms and mean images computed' + pol_txt)

            ## check if variation of rms/mean maps is significant enough:
            #       check_rmsmap() sets img.use_rms_map
            #       check_meanmap() sets img.mean_map_type
            if pol == 'I':
                if opts.rms_map is None and img.use_rms_map is None:
                    if do_adapt and len(isl_pos) > 0:
                        # Always use 2d map if there is at least one bright
                        # source and adaptive scaling is desired
                        img.use_rms_map = True
                    else:
                        self.check_rmsmap(img, rms)
                elif opts.rms_map != None:
                    img.use_rms_map = opts.rms_map
                if img.use_rms_map is False:
                    mylogger.userinfo(mylog, 'Using constant background rms')
                else:
                    mylogger.userinfo(mylog, 'Using 2D map for background rms')

                if opts.mean_map == 'default' and img.mean_map_type is None:
                    self.check_meanmap(img, rms)
                elif opts.mean_map != 'default':
                    img.mean_map_type = opts.mean_map
                if img.mean_map_type != 'map':
                    mylogger.userinfo(mylog, 'Using constant background mean')
                else:
                    mylogger.userinfo(mylog, 'Using 2D map for background mean')

          ## if rms map is insignificant, or rms_map==False use const value
          if img.use_rms_map is False:
            if opts.rms_value == None:
              rms[:]  = crmss[ipol]
            else:
              rms[:]  = opts.rms_value
            mylogger.userinfo(mylog, 'Value of background rms' + pol_txt,
                              '%.5f Jy/beam' % rms[0][0])
          else:
            rms_min = N.nanmin(rms)
            rms_max = N.nanmax(rms)
            mylogger.userinfo(mylog, 'Min/max values of background rms map' + pol_txt,
                              '(%.5f, %.5f) Jy/beam' % (rms_min, rms_max))

          if img.mean_map_type != 'map':
            if opts.mean_map == 'zero':
                val = 0.0
            else:
                val = img.clipped_mean
            mean[:] = val
            mylogger.userinfo(mylog, 'Value of background mean' + pol_txt,
                              str(round(val,5))+' Jy/beam')
          else:
            mean_min = N.nanmin(mean)
            mean_max = N.nanmax(mean)
            mylogger.userinfo(mylog, 'Min/max values of background mean map' + pol_txt,
                              '(%.5f, %.5f) Jy/beam' % (mean_min, mean_max))

          if pol == 'I':
            # Apply mask to mean_map and rms_map by setting masked values to NaN
            if isinstance(mask, N.ndarray):
                pix_masked = N.where(mask == True)
                mean[pix_masked] = N.nan
                rms[pix_masked] = N.nan
            img.mean = mean; img.rms  = rms
            if opts.savefits_rmsim or opts.output_all:
              if img.waveletimage:
                  resdir = img.basedir + '/wavelet/background/'
              else:
                  resdir = img.basedir + '/background/'
              if not os.path.exists(resdir): os.makedirs(resdir)
              func.write_image_to_file(img.use_io, img.imagename + '.rmsd_I.fits', rms, img, resdir)
              mylog.info('%s %s' % ('Writing ', resdir+img.imagename+'.rmsd_I.fits'))
            if opts.savefits_meanim or opts.output_all:
              if img.waveletimage:
                  resdir = img.basedir + '/wavelet/background/'
              else:
                  resdir = img.basedir + '/background/'
              if not os.path.exists(resdir): os.makedirs(resdir)
              func.write_image_to_file(img.use_io, img.imagename + '.mean_I.fits', mean, img, resdir)
              mylog.info('%s %s' % ('Writing ', resdir+img.imagename+'.mean_I.fits'))
            if opts.savefits_normim or opts.output_all:
              if img.waveletimage:
                  resdir = img.basedir + '/wavelet/background/'
              else:
                  resdir = img.basedir + '/background/'
              if not os.path.exists(resdir): os.makedirs(resdir)
              zero_pixels = N.where(rms <= 0.0)
              rms_nonzero = rms.copy()
              rms_nonzero[zero_pixels] = N.NaN
              func.write_image_to_file(img.use_io, img.imagename + '.norm_I.fits', (img.ch0-mean)/rms_nonzero, img, resdir)
              mylog.info('%s %s' % ('Writing ', resdir+img.imagename+'.norm_I.fits'))
          else:
            img.mean_QUV.append(mean); img.rms_QUV.append(rms)

        img.completed_Ops.append('rmsimage')
        return img

    def check_rmsmap(self, img, rms):
        """Calculates the statistics of the rms map and decides, when
        rms_map=None, whether to take the map (if variance
        is significant) or a constant value
        """
    	from math import sqrt

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Rmsimage.Checkrms  ")
        cdelt = img.wcs_obj.acdelt[:2]
    	bm = (img.beam[0], img.beam[1])
    	fw_pix = sqrt(N.product(bm)/abs(N.product(cdelt)))
    	if img.masked:
    	    unmasked = N.where(~img.mask)
            stdsub = N.std(rms[unmasked])
            maxrms = N.max(rms[unmasked])
        else:
            stdsub = N.std(rms)
            maxrms = N.max(rms)

    	rms_expect = img.clipped_rms/sqrt(2)/img.rms_box[0]*fw_pix
        mylog.debug('%s %10.6f %s' % ('Standard deviation of rms image = ', stdsub*1000.0, 'mJy'))
        mylog.debug('%s %10.6f %s' % ('Expected standard deviation = ', rms_expect*1000.0, 'mJy'))
    	if stdsub > 1.1*rms_expect:
            img.use_rms_map = True
            mylogger.userinfo(mylog, 'Variation in rms image significant')
        else:
            img.use_rms_map = False
            mylogger.userinfo(mylog, 'Variation in rms image not significant')

        return img

    def check_meanmap(self, img, mean):
        """Calculates the statistics of the mean map and decides, when
        mean_map=None, whether to take the map (if variance
        is significant) or a constant value
        """
    	from math import sqrt

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Rmsimage.Checkmean ")
        cdelt = img.wcs_obj.acdelt[:2]
        bm = (img.beam[0], img.beam[1])
        fw_pix = sqrt(N.product(bm)/abs(N.product(cdelt)))
    	if img.masked:
            unmasked = N.where(~img.mask)
            stdsub = N.std(mean[unmasked])
            maxmean = N.max(mean[unmasked])
        else:
            stdsub = N.std(mean)
            maxmean = N.max(mean)
        rms_expect = img.clipped_rms/img.rms_box[0]*fw_pix
        mylog.debug('%s %10.6f %s' % ('Standard deviation of mean image = ', stdsub*1000.0, 'mJy'))
        mylog.debug('%s %10.6f %s' % ('Expected standard deviation = ', rms_expect*1000.0, 'mJy'))

        # For mean map, use a higher threshold than for the rms map, as radio images
        # should rarely, if ever, have significant variations in the mean
        if stdsub > 5.0*rms_expect:
          img.mean_map_type = 'map'
          mylogger.userinfo(mylog, 'Variation in mean image significant')
        else:
          if img.confused:
            img.mean_map_type = 'zero'
          else:
            img.mean_map_type = 'const'
          mylogger.userinfo(mylog, 'Variation in mean image not significant')

        return img


    def calculate_maps(self, img, data, mean, rms, mask, map_opts, do_adapt,
                       bright_pt_coords=[], rms_box2=None,
                       logname=None, ncores=None):
        """Calls map_2d and checks for problems"""
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Rmsimage.Calcmaps ")
        rms_ok = False
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Rmsimage.Calcmaps ")
        opts = img.opts
        while not rms_ok:
            self.map_2d(data, mean, rms, mask, *map_opts, do_adapt=do_adapt,
                        bright_pt_coords=bright_pt_coords, rms_box2=rms_box2,
                        logname=logname, ncores=ncores)
            if N.any(rms < 0.0):
                rms_ok = False
                if (opts.rms_box_bright is None and do_adapt) or (opts.rms_box is None and not do_adapt):
                    # Increase box by 20%
                    if do_adapt:
                        new_width = int(img.rms_box_bright[0]*1.2)
                        if new_width == img.rms_box_bright[0]:
                            new_width = img.rms_box_bright[0] + 1
                        new_step = int(new_width/3.0)
                        img.rms_box_bright = (new_width, new_step)
                        if img.rms_box_bright[0] > min(img.ch0.shape)/4.0:
                            mylogger.userinfo(mylog, 'Size of rms_box_bright larger than 1/4 of image size')
                            mylogger.userinfo(mylog, 'Using constant background rms and mean')
                            img.use_rms_map = False
                            img.rms_box = img.rms_box_bright
                            img.mean_map_type = 'const'
                            rms_ok = True
                        else:
                            map_opts = (opts.kappa_clip, img.rms_box_bright, opts.spline_rank)
                    else:
                        new_width = int(img.rms_box[0]*1.2)
                        if new_width == img.rms_box[0]:
                            new_width = img.rms_box[0] + 1
                        new_step = int(new_width/3.0)
                        img.rms_box = (new_width, new_step)
                        if img.rms_box[0] > min(img.ch0.shape)/4.0:
                            mylogger.userinfo(mylog, 'Size of rms_box larger than 1/4 of image size')
                            mylogger.userinfo(mylog, 'Using constant background rms and mean')
                            img.use_rms_map = False
                            img.mean_map_type = 'const'
                            rms_ok = True
                        else:
                            map_opts = (opts.kappa_clip, img.rms_box, opts.spline_rank)

                else:
                    # User has specified box size, use order=1 to prevent negatives
                    if opts.spline_rank > 1:
                        mylog.warning('Negative values found in rms map interpolated with spline_rank = %i' % opts.spline_rank)
                        mylog.warning('Using spline_rank = 1 (bilinear interpolation) instead')
                        if do_adapt:
                            map_opts = (opts.kappa_clip, img.rms_box_bright, 1)
                        else:
                            map_opts = (opts.kappa_clip, img.rms_box, 1)
            else:
                rms_ok = True

        return mean, rms


    def map_2d(self, arr, out_mean, out_rms, mask=False,
               kappa=3, box=None, interp=1, do_adapt=False,
               bright_pt_coords=None, rms_box2=None, logname='', ncores=None):
        """Calculate mean&rms maps and store them into provided arrays

        Parameters:
        arr: 2D array with data
        out_mean, out_rms: 2D arrays where to store calculated maps
        mask: mask
        kappa: clipping value for rms/mean calculations
        box: tuple of (box_size, box_step) for calculating map
        rms_box2 = large-scale box size
        interp: order of interpolating spline used to interpolate
                calculated map
        do_adapt: use adaptive binning
        """
        mask_small = mask
        axes, mean_map1, rms_map1 = self.rms_mean_map(arr, mask_small, kappa, box, ncores)
        ax = map(self.remap_axis, arr.shape, axes)
        ax = N.meshgrid(*ax[-1::-1])
        pt_src_scale = box[0]
        if do_adapt:
            out_rms2 = N.zeros(rms_map1.shape)
            out_mean2 = N.zeros(rms_map1.shape)
            # Generate rms/mean maps on large scale
            box2 = rms_box2
            axes2, mean_map2, rms_map2 = self.rms_mean_map(arr, mask, kappa, box2, ncores)

            # Interpolate to get maps on small scale grid
            axes2mod = axes2[:]
            axes2mod[0] = axes2[0]/arr.shape[0]*mean_map1.shape[0]
            axes2mod[1] = axes2[1]/arr.shape[1]*mean_map1.shape[1]
            ax2 = map(self.remap_axis, out_rms2.shape, axes2mod)
            ax2 = N.meshgrid(*ax2[-1::-1])
            nd.map_coordinates(rms_map2,  ax2[-1::-1], order=interp, output=out_rms2)
            nd.map_coordinates(mean_map2, ax2[-1::-1], order=interp, output=out_mean2)
            rms_map = out_rms2
            mean_map = out_mean2

            # For each bright source, find nearest points and weight them towards
            # the small scale maps.
            xscale = float(arr.shape[0])/float(out_rms2.shape[0])
            yscale = float(arr.shape[1])/float(out_rms2.shape[1])
            scale = [xscale, yscale]
            size = 15
            for bright_pt in bright_pt_coords:
                bbox, src_center = self.make_bright_src_bbox(bright_pt, scale, size, out_rms2.shape)
                bbox_xsize = bbox[0].stop-bbox[0].start
                bbox_ysize = bbox[1].stop-bbox[1].start
                src_center[0] -= bbox[0].start
                src_center[1] -= bbox[1].start
                weights = N.ones((bbox_xsize, bbox_ysize))

                # Taper weights to zero where small-scale value is within a factor of
                # 2 of large-scale value. Use distance to center of the box
                # to determine taper value. This tapering prevents the use of the
                # small-scale box beyond the range of artifacts.
                low_vals_ind = N.where(rms_map1[bbox]/out_rms2[bbox] < 2.0)
                if len(low_vals_ind[0]) > 0:
                    dist_to_cen = []
                    for (x,y) in zip(low_vals_ind[0],low_vals_ind[1]):
                        dist_to_cen.append(N.sqrt( (x-src_center[0])**2 +
                                           (y-src_center[1])**2 ))
                    med_dist_to_cen = N.min(dist_to_cen)
                    for x in range(bbox_xsize):
                        for y in range(bbox_ysize):
                            dist_to_cen = N.sqrt( (x-src_center[0])**2 +
                                               (y-src_center[1])**2 )
                            if dist_to_cen >= med_dist_to_cen:
                                weights[x,y] = 1.0 - dist_to_cen/N.sqrt(bbox_xsize**2+bbox_ysize**2)*2.0
                rms_map[bbox] = rms_map1[bbox]*weights + out_rms2[bbox]*(1.0-weights)
                mean_map[bbox] = mean_map1[bbox]*weights + out_mean2[bbox]*(1.0-weights)
        else:
            rms_map = rms_map1
            mean_map = mean_map1

        # Interpolate to image coords
        # Check whether default (cubic) or user-specified order causes problems.
        # If so, use order=1.
        mylog = mylogger.logging.getLogger(logname+"Rmsimage")
        nd.map_coordinates(rms_map,  ax[-1::-1], order=interp, output=out_rms)
        nd.map_coordinates(mean_map, ax[-1::-1], order=interp, output=out_mean)

        # Apply mask to mean_map and rms_map by setting masked values to NaN
        if isinstance(mask, N.ndarray):
            pix_masked = N.where(mask == True)
            out_mean[pix_masked] = N.nan
            out_rms[pix_masked] = N.nan

    def rms_mean_map(self, arr, mask=False, kappa=3, box=None, ncores=None):
        """Calculate map of the mean/rms values

        Parameters:
        arr:  2D array with data
        mask: mask
        kappa: clipping for calculating rms/mean within each box
        box: box parameters (box_size, box_step)

        Returns:
        axes: list of 2 arrays with coordinates of boxes alongside each axis
        mean_map: map of mean values
        rms_map: map of rms values

        Description:
        This function calculates clipped mean and rms maps for the array.
        The algorithm is a moving-window algorithm, where mean&rms are
        calculated within a window of a size (box_size * box_size), and the
        window is stepped withing the image by steps of box_steps.

        Special care is taken for the borders of the image -- outer borders
        (where box doesn't fit properly) are given one extra round with a box
        applied to the border of the image. Additionally outer values are
        extrapolated to cover whole image size, to simplify further processing.

        See also routine 'remap_axes' for 'inverting' axes array

        Example:
        for an input image of 100x100 pixels calling rms_mean_map with default
        box parameters (50, 25) will result in the following:

        axes = [array([  0. ,  24.5,  49.5,  74.5,  99. ]),
                array([  0. ,  24.5,  49.5,  74.5,  99. ])]

        mean_map = <5x5 array>
        rms_map  = <5x5 array>

        rms_map[1,1] is calculated for  arr[0:50, 0:50]
        rms_map[2,1] is calculated for  arr[25:75, 0:50]
        ...etc...
        rms_map[0,0] is extrapolated as .5*(rms_map[0,1] + rms_map[1,0])
        rms_map[0,1] is extrapolated as rms_map[1,1]
        """
        mylog = mylogger.logging.getLogger("PyBDSM.RmsMean")
        if box is None:
            box = (50, 25)
        if box[0] < box[1]:
            raise RuntimeError('Box size is less than step size.')

        # Some math first: boxcount is number of boxes alongsize each axis,
        # bounds is non-zero for axes which have extra pixels beyond last box
        BS, SS = box
        imgshape = N.array(arr.shape)

        # If boxize is less than 10% of image, use simple extrapolation to
        # derive the edges of the mean and rms maps; otherwise, use padded
        # versions of arr and mask to derive the mean and rms maps
        if float(BS)/float(imgshape[0]) < 0.1 and \
                float(BS)/float(imgshape[1]) < 0.1:
            use_extrapolation = True
        else:
            use_extrapolation = False

        if use_extrapolation:
            boxcount = 1 + (imgshape - BS)/SS
            bounds   = N.asarray((boxcount-1)*SS + BS < imgshape, dtype=int)
            mapshape = 2 + boxcount + bounds
        else:
            boxcount = 1 + imgshape/SS
            bounds   = N.asarray((boxcount-1)*SS < imgshape, dtype=int)
            mapshape = boxcount + bounds
            pad_border_size = int(BS/2.0)
            new_shape = (arr.shape[0] + 2*pad_border_size, arr.shape[1]
                         + 2*pad_border_size)
            arr_pad = self.pad_array(arr, new_shape)
            if mask == None:
                mask_pad = None
            else:
                mask_pad = self.pad_array(mask, new_shape)

        # Make arrays for calculated data
        mean_map = N.zeros(mapshape, dtype=float)
        rms_map  = N.zeros(mapshape, dtype=float)
        axes     = [N.zeros(len, dtype=float) for len in mapshape]

        # Step 1: internal area of the image
        # Make a list of coordinates to send to process_mean_rms_maps()
        coord_list = []
        ind_list = []
        for i in range(boxcount[0]):
            for j in range(boxcount[1]):
                if use_extrapolation:
                    coord_list.append((i+1, j+1))
                else:
                    coord_list.append((i, j))
                ind_list.append([i*SS, i*SS+BS, j*SS, j*SS+BS])

        # Now call the parallel mapping function. Returns a list of [mean, rms]
        # for each coordinate.
        if use_extrapolation:
            cm_cr_list = mp.parallel_map(func.eval_func_tuple,
                    itertools.izip(itertools.repeat(self.process_mean_rms_maps),
                    ind_list, itertools.repeat(mask), itertools.repeat(arr),
                    itertools.repeat(kappa)), numcores=ncores)
        else:
            cm_cr_list = mp.parallel_map(func.eval_func_tuple,
                    itertools.izip(itertools.repeat(self.process_mean_rms_maps),
                    ind_list, itertools.repeat(mask_pad), itertools.repeat(arr_pad),
                    itertools.repeat(kappa)), numcores=ncores)

        for i, co in enumerate(coord_list):
            cm, cr = cm_cr_list[i]
            mean_map[co] = cm
            rms_map[co] = cr

        # Check if all regions have too few unmasked pixels
        if mask != None and N.size(N.where(mean_map != N.inf)) == 0:
            raise RuntimeError("No unmasked regions from which to determine "\
                         "mean and rms maps")

        # Step 2: borders of the image
        if bounds[0]:
            coord_list = []
            ind_list = []
            for j in range(boxcount[1]):
                if use_extrapolation:
                    coord_list.append((-2, j+1))
                    ind_list.append([-BS, arr.shape[0], j*SS,j*SS+BS])
                else:
                    coord_list.append((-1, j))
                    ind_list.append([-BS, arr_pad.shape[0], j*SS,j*SS+BS])
            if use_extrapolation:
                cm_cr_list = mp.parallel_map(func.eval_func_tuple,
                        itertools.izip(itertools.repeat(self.process_mean_rms_maps),
                        ind_list, itertools.repeat(mask), itertools.repeat(arr),
                        itertools.repeat(kappa)), numcores=ncores)
            else:
                cm_cr_list = mp.parallel_map(func.eval_func_tuple,
                        itertools.izip(itertools.repeat(self.process_mean_rms_maps),
                        ind_list, itertools.repeat(mask_pad), itertools.repeat(arr_pad),
                        itertools.repeat(kappa)), numcores=ncores)

            for i, co in enumerate(coord_list):
                cm, cr = cm_cr_list[i]
                mean_map[co] = cm
                rms_map[co] = cr


        if bounds[1]:
            coord_list = []
            ind_list = []
            for i in range(boxcount[0]):
                if use_extrapolation:
                    coord_list.append((i+1, -2))
                    ind_list.append([i*SS,i*SS+BS, -BS,arr.shape[1]])
                else:
                    coord_list.append((i, -1))
                    ind_list.append([i*SS,i*SS+BS, -BS,arr_pad.shape[1]])
            if use_extrapolation:
                cm_cr_list = mp.parallel_map(func.eval_func_tuple,
                        itertools.izip(itertools.repeat(self.process_mean_rms_maps),
                        ind_list, itertools.repeat(mask), itertools.repeat(arr),
                        itertools.repeat(kappa)), numcores=ncores)
            else:
                cm_cr_list = mp.parallel_map(func.eval_func_tuple,
                        itertools.izip(itertools.repeat(self.process_mean_rms_maps),
                        ind_list, itertools.repeat(mask_pad), itertools.repeat(arr_pad),
                        itertools.repeat(kappa)), numcores=ncores)

            for i, co in enumerate(coord_list):
                cm, cr = cm_cr_list[i]
                mean_map[co] = cm
                rms_map[co] = cr

        if bounds.all():
                if use_extrapolation:
                    ind = [-BS,arr.shape[0], -BS,arr.shape[1]]
                    self.for_masked(mean_map, rms_map, mask, arr, ind,
                                    kappa, [-2, -2])
                else:
                    ind = [-BS,arr_pad.shape[0], -BS,arr_pad.shape[1]]
                    self.for_masked(mean_map, rms_map, mask_pad, arr_pad, ind,
                                    kappa, [-1, -1])

        # Step 3: correct(extrapolate) borders of the image
        def correct_borders(map):
            map[0, :] = map[1, :]
            map[:, 0] = map[:, 1]
            map[-1, :] = map[-2, :]
            map[:, -1] = map[:, -2]

            map[0,0] = (map[1,0] + map[0, 1])/2.
            map[-1,0] = (map[-2, 0] + map[-1, 1])/2.
            map[0, -1] = (map[0, -2] + map[1, -1])/2.
            map[-1,-1] = (map[-2, -1] + map[-1, -2])/2.

        if use_extrapolation:
            correct_borders(mean_map)
            correct_borders(rms_map)

        # Step 4: fill in coordinate axes
        for i in range(2):
            if use_extrapolation:
                axes[i][1:boxcount[i]+1] = (N.arange(boxcount[i])*SS
                                            + BS/2. - .5)
                if bounds[i]:
                    axes[i][-2] = imgshape[i] - BS/2. - .5
            else:
                axes[i][0:boxcount[i]] = N.arange(boxcount[i])*SS - .5
                if bounds[i]:
                    axes[i][-2] = imgshape[i] - .5
            axes[i][-1] = imgshape[i] - 1

        # Step 5: fill in boxes with < 5 unmasked pixels (set to values of
        # N.inf)
        unmasked_boxes = N.where(mean_map != N.inf)
        if N.size(unmasked_boxes,1) < mapshape[0]*mapshape[1]:
            mean_map = self.fill_masked_regions(mean_map)
            rms_map = self.fill_masked_regions(rms_map)

        return axes, mean_map, rms_map


    def process_mean_rms_maps(self, ind, mask, arr, kappa):
        """Finds mean and rms for one region of an input arr"""
        cm, cr = self.for_masked_mp(mask, arr, ind,
                        kappa)
        return cm, cr


    def fill_masked_regions(self, themap, magic=N.inf):
        """Fill masked regions (defined where values == magic) in themap.
        """
        masked_boxes = N.where(themap == magic) # locations of masked regions
        for i in range(N.size(masked_boxes,1)):
            num_unmasked = 0
            x, y = masked_boxes[0][i], masked_boxes[1][i]
            delx = dely = 1
            while num_unmasked == 0:
                x1 = x - delx
                if x1 < 0: x1 = 0
                x2 = x + 1 + delx
                if x2 > themap.shape[0]: x2 = themap.shape[0]
                y1 = y - dely
                if y1 < 0: y1 = 0
                y2 = y + 1 + dely
                if y2 > themap.shape[1]: y2 = themap.shape[1]

                cutout = themap[x1:x2, y1:y2].ravel()
                goodcutout = cutout[cutout != magic]
                num_unmasked = N.alen(goodcutout)
                if num_unmasked > 0:
                    themap[x, y] = N.nansum(goodcutout)/float(len(goodcutout))
                delx += 1
                dely += 1
        themap[N.where(N.isnan(themap))] = 0.0
        return themap

    def pad_array(self, arr, new_shape):
        """Returns a padded array by mirroring around the edges."""
        # Assume that padding is the same for both axes and is equal
        # around all edges.
        half_size = (new_shape[0] - arr.shape[0]) / 2
        arr_pad = N.zeros( (new_shape), dtype=arr.dtype )

        # left band
        band = arr[:half_size, :]
        arr_pad[:half_size, half_size:-half_size] =  N.flipud( band )

        # right band
        band = arr[-half_size:, :]
        arr_pad[-half_size:, half_size:-half_size] = N.flipud( band )

        # bottom band
        band = arr[:, :half_size]
        arr_pad[half_size:-half_size, :half_size] = N.fliplr( band )

        # top band
        band = arr[:, -half_size:]
        arr_pad[half_size:-half_size, -half_size:] =  N.fliplr( band )

        # central band
        arr_pad[half_size:-half_size, half_size:-half_size] = arr

        # bottom left corner
        band = arr[:half_size,:half_size]
        arr_pad[:half_size,:half_size] = N.flipud(N.fliplr(band))

        # top right corner
        band = arr[-half_size:,-half_size:]
        arr_pad[-half_size:,-half_size:] = N.flipud(N.fliplr(band))

        # top left corner
        band = arr[:half_size,-half_size:]
        arr_pad[:half_size,-half_size:] = N.flipud(N.fliplr(band))

        # bottom right corner
        band = arr[-half_size:,:half_size]
        arr_pad[-half_size:,:half_size] = N.flipud(N.fliplr(band))

        return arr_pad

    def for_masked(self, mean_map, rms_map, mask, arr, ind, kappa, co):

        bstat = _cbdsm.bstat
        a, b, c, d = ind; i, j = co
        if mask == None:
          m, r, cm, cr, cnt = bstat(arr[a:b, c:d], mask, kappa)
          if cnt > 198: cm = m; cr = r
          mean_map[i, j], rms_map[i, j] = cm, cr
        else:
          pix_unmasked = N.where(mask[a:b, c:d] == False)
          npix_unmasked = N.size(pix_unmasked,1)
          if npix_unmasked > 20: # find clipped mean/rms
            m, r, cm, cr, cnt = bstat(arr[a:b, c:d], mask[a:b, c:d], kappa)
            if cnt > 198: cm = m; cr = r
            mean_map[i, j], rms_map[i, j] = cm, cr
          else:
            if npix_unmasked > 5: # just find simple mean/rms
              cm = N.mean(arr[pix_unmasked])
              cr = N.std(arr[pix_unmasked])
              mean_map[i, j], rms_map[i, j] = cm, cr
            else: # too few unmasked pixels --> set mean/rms to inf
              mean_map[i, j], rms_map[i, j] = N.inf, N.inf


    def for_masked_mp(self, mask, arr, ind, kappa):

        bstat = _cbdsm.bstat
        a, b, c, d = ind
        if mask == None:
          m, r, cm, cr, cnt = bstat(arr[a:b, c:d], mask, kappa)
          if cnt > 198: cm = m; cr = r
        else:
          pix_unmasked = N.where(mask[a:b, c:d] == False)
          npix_unmasked = N.size(pix_unmasked,1)
          if npix_unmasked > 20: # find clipped mean/rms
            m, r, cm, cr, cnt = bstat(arr[a:b, c:d], mask[a:b, c:d], kappa)
            if cnt > 198: cm = m; cr = r
          else:
            if npix_unmasked > 5: # just find simple mean/rms
              cm = N.mean(arr[pix_unmasked])
              cr = N.std(arr[pix_unmasked])
            else: # too few unmasked pixels --> set mean/rms to inf
              cm = N.inf
              cr = N.inf

        return cm, cr


    def remap_axis(self, size, arr):
        """Invert axis mapping done by rms_mean_map

        rms_mean_map 'compresses' axes by returning short arrays with
        coordinades of the boxes. This routine 'inverts' this compression
        by calculating coordinates of each pixel of the original array
        within compressed array.

        Parameters:
        size: size of the original (and resulting) array
        arr : 'compressed' axis array from rms_mean_map

        Example:
        the following 'compressed' axis (see example in rms_mean_map):

           ax = array([  0. ,  24.5,  49.5,  74.5,  99. ])

        will be remapped as:

        print remap_axis(100, ax)
        [ 0.          0.04081633  0.08163265  0.12244898 ....
           ...............................................
          3.91836735  3.95918367  4.        ]

        which means that pixel 0 in the original image corresponds to pixels
        0 in the rms/mean_map array (which is 5x5 array).
        pixel 1 of the original image has coordinate of 0.04081633 in the
        compressed image (e.g. it has no exact counterpart, and it's value
        should be obtained by interpolation)
        """
        from math import floor, ceil
        res = N.zeros(size, dtype=float)

        for i in range(len(arr) - 1):
            i1 = arr[i]
            i2 = arr[i+1]
            t = N.arange(ceil(i1), floor(i2)+1, dtype=float)
            res[ceil(i1):floor(i2)+1] = i + (t-i1)/(i2-i1)

        return res

    def make_bright_src_bbox(self, coord, scale, size, shape):
        """Returns bbox given coordinates of center and scale"""
        xindx = int(coord[0]/scale[0])
        yindx = int(coord[1]/scale[1])
        xlow = xindx - int(size/2.0)
        if xlow < 0:
            xlow = 0
        xhigh = xindx + int(size/2.0) + 1
        if xhigh > shape[0]:
            xhigh = shape[0]
        ylow = yindx - int(size/2.0)
        if ylow < 0:
            ylow = 0
        yhigh = yindx + int(size/2.0) + 1
        if yhigh > shape[1]:
            yhigh = shape[1]

        src_center = [xindx, yindx]
        return [slice(xlow, xhigh, None), slice(ylow, yhigh, None)], src_center

    def output_rmsbox_size(self, img):
        """Prints rms/mean box size"""
        opts = img.opts
        do_adapt = opts.adaptive_rms_box
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"RMSimage")
        if (opts.rms_map is not False) or (opts.mean_map not in ['zero', 'const']):
          if do_adapt:
              if opts.rms_box_bright is None:
                  mylogger.userinfo(mylog, 'Derived rms_box (box size, step size)',
                                    '(' + str(img.rms_box_bright[0]) + ', ' +
                                    str(img.rms_box_bright[1]) + ') pixels (small scale)')
              else:
                  mylogger.userinfo(mylog, 'Using user-specified rms_box',
                                    '(' + str(img.rms_box_bright[0]) + ', ' +
                                    str(img.rms_box_bright[1]) + ') pixels (small scale)')
              if opts.rms_box is None:
                  mylogger.userinfo(mylog, 'Derived rms_box (box size, step size)',
                                    '(' + str(img.rms_box[0]) + ', ' +
                                    str(img.rms_box[1]) + ') pixels (large scale)')
              else:
                  mylogger.userinfo(mylog, 'Using user-specified rms_box',
                                    '(' + str(img.rms_box[0]) + ', ' +
                                    str(img.rms_box[1]) + ') pixels (large scale)')
          else:
              if opts.rms_box is None:
                  mylogger.userinfo(mylog, 'Derived rms_box (box size, step size)',
                                '(' + str(img.rms_box[0]) + ', ' +
                                str(img.rms_box[1]) + ') pixels')
              else:
                  mylogger.userinfo(mylog, 'Using user-specified rms_box',
                                    '(' + str(img.rms_box[0]) + ', ' +
                                    str(img.rms_box[1]) + ') pixels')
