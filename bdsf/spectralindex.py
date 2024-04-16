"""Module Spectral index.

   This module calculates spectral indices for Gaussians and sources for a multichannel cube.

"""
from __future__ import print_function
from __future__ import absolute_import

import numpy as N
from .image import Op
from . import mylogger
from copy import deepcopy as cp
from . import functions as func
from . import statusbar


class Op_spectralindex(Op):
    """Computes spectral index of every gaussian and every source.

    First do a quick fit to all channels to determine whether averaging over
    frequency is needed to obtain desired SNR (set by img.opts.specind_snr).
    This averaging should be done separately for both Gaussians and
    sources. For S and C sources, averaging only needs to be done once
    (as the sources have only one Gaussian).

    For M sources, averaging is needed twice: once to obtain the desired
    SNR for the faintest Gaussian in the source, and once to obtain the
    desired SNR for the source as a whole.

    If averaging is needed for a given source, don't let the
    number of resulting channels fall below 2. If it is not possible
    to obtain the desired SNR in 2 or more channels, set spec_indx of
    Gaussian/source to NaN.

   """

    def __call__(self, img):
        global bar1

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"SpectIndex")
        img.mylog = mylog
        if img.opts.spectralindex_do:
            mylogger.userinfo(mylog, '\nExtracting spectral indices for all ch0 sources')
            shp = img.image_arr.shape
            if shp[1] > 1:
                # calc freq, beam_spectrum for nchan channels
                self.freq_beamsp_unav(img)
                sbeam = img.beam_spectrum
                freqin = img.freq

                # calc initial channel flags if needed
                iniflags = self.iniflag(img)
                img.specind_iniflags = iniflags
                good_chans = N.where(iniflags == False)
                unav_image = img.image_arr[0][good_chans]
                unav_freqs = freqin[good_chans]
                nmax_to_avg = img.opts.specind_maxchan
                nchan = unav_image.shape[0]
                mylogger.userinfo(mylog, 'Number of channels remaining after initial flagging', str(nchan))
                if nmax_to_avg == 0:
                    nmax_to_avg = nchan

                # calculate the rms map of each unflagged channel
                bar1 = statusbar.StatusBar('Determing rms for channels in image ..... : ', 0, nchan)
                if img.opts.quiet == False:
                    bar1.start()
                rms_spec = self.rms_spectrum(img, unav_image) # bar1 updated here

                bar2 = statusbar.StatusBar('Calculating spectral indices for sources  : ', 0, img.nsrc)
                c_wts = img.opts.collapse_wt
                snr_desired = img.opts.specind_snr

                if img.opts.quiet == False and img.opts.verbose_fitting == False:
                    bar2.start()
                for src in img.sources:
                    isl = img.islands[src.island_id]
                    isl_bbox = isl.bbox

                    # Fit each channel with ch0 Gaussian(s) of the source,
                    # allowing only the normalization to vary.
                    chan_images = unav_image[:, isl_bbox[0], isl_bbox[1]]
                    chan_rms = rms_spec[:, isl_bbox[0], isl_bbox[1]]
                    beamlist = img.beam_spectrum
                    unavg_total_flux, e_unavg_total_flux = self.fit_channels(img, chan_images, chan_rms, src, beamlist)

                    # Check for upper limits and mask. gaus_mask is array of (N_channels x N_gaussians)
                    # and is True if measured flux is upper limit. n_good_chan_per_gaus is array of N_gaussians
                    # that gives number of unmasked channels for each Gaussian.
                    gaus_mask, n_good_chan_per_gaus = self.mask_upper_limits(unavg_total_flux, e_unavg_total_flux, snr_desired)

                    # Average if needed and fit again
                    # First find flux of faintest Gaussian of source and use it to estimate rms_desired
                    gflux = []
                    for g in src.gaussians:
                        gflux.append(g.peak_flux)
                    rms_desired = min(gflux)/snr_desired
                    total_flux = unavg_total_flux
                    e_total_flux = e_unavg_total_flux
                    freq_av = unav_freqs
                    nchan = chan_images.shape[0]
                    nchan_prev = nchan
                    while min(n_good_chan_per_gaus) < 2 and nchan > 2:
                        avimages, beamlist, freq_av, crms_av = self.windowaverage_cube(chan_images, rms_desired, chan_rms,
                                    c_wts, sbeam, freqin, nmax_to_avg=nmax_to_avg)
                        total_flux, e_total_flux = self.fit_channels(img, avimages, crms_av, src, beamlist)
                        gaus_mask, n_good_chan_per_gaus = self.mask_upper_limits(total_flux, e_total_flux, snr_desired)
                        nchan = avimages.shape[0]
                        if nchan == nchan_prev:
                            break
                        nchan_prev = nchan
                        rms_desired *= 0.8

                    # Now fit Gaussian fluxes to obtain spectral indices.
                    # Only fit if there are detections (at specified sigma threshold)
                    # in at least two bands. If not, don't fit and set spec_indx
                    # and error to NaN.
                    for ig, gaussian in enumerate(src.gaussians):
                        npos = len(N.where(total_flux[:, ig] > 0.0)[0])
                        if img.opts.verbose_fitting:
                            if img.opts.flagchan_snr:
                                print('Gaussian #%i : averaged to %i channels, of which %i meet SNR criterion' % (gaussian.gaus_num,
                                       len(total_flux[:, ig]), n_good_chan_per_gaus[ig]))
                            else:
                                print('Gaussian #%i : averaged to %i channels, all of which will be used' % (gaussian.gaus_num,
                                       len(total_flux[:, ig])))
                        if (img.opts.flagchan_snr and n_good_chan_per_gaus[ig] < 2) or npos < 2:
                            gaussian.spec_indx = N.NaN
                            gaussian.e_spec_indx = N.NaN
                            gaussian.spec_norm = N.NaN
                            gaussian.specin_flux = [N.NaN]
                            gaussian.specin_fluxE = [N.NaN]
                            gaussian.specin_freq = [N.NaN]
                            gaussian.specin_freq0 = N.NaN
                        else:
                            if img.opts.flagchan_snr:
                                good_fluxes_ind = N.where(gaus_mask[:, ig] == False)
                            else:
                                good_fluxes_ind = range(len(freq_av))
                            fluxes_to_fit = total_flux[:, ig][good_fluxes_ind]
                            e_fluxes_to_fit = e_total_flux[:, ig][good_fluxes_ind]
                            freqs_to_fit = freq_av[good_fluxes_ind]

                            fit_res = self.fit_specindex(freqs_to_fit, fluxes_to_fit, e_fluxes_to_fit)
                            gaussian.spec_norm, gaussian.spec_indx, gaussian.e_spec_indx = fit_res
                            gaussian.specin_flux = fluxes_to_fit.tolist()
                            gaussian.specin_fluxE = e_fluxes_to_fit.tolist()
                            gaussian.specin_freq = freqs_to_fit.tolist()
                            gaussian.specin_freq0 = N.median(freqs_to_fit)

                    # Next fit total source fluxes for spectral index.
                    if len(src.gaussians) > 1:
                        # First, check unaveraged SNRs for total source.
                        src_total_flux = N.zeros((chan_images.shape[0], 1))
                        src_e_total_flux = N.zeros((chan_images.shape[0], 1))
                        src_total_flux[:,0] = N.sum(unavg_total_flux, 1) # sum over all Gaussians in source to obtain total fluxes in each channel
                        src_e_total_flux[:,0] = N.sqrt(N.sum(N.power(e_unavg_total_flux, 2.0), 1))
                        src_mask, n_good_chan = self.mask_upper_limits(src_total_flux, src_e_total_flux, snr_desired)

                        # Average if needed and fit again
                        rms_desired = src.peak_flux_max/snr_desired
                        total_flux = unavg_total_flux
                        e_total_flux = e_unavg_total_flux
                        freq_av = unav_freqs
                        nchan = chan_images.shape[0]
                        nchan_prev = nchan
                        while n_good_chan < 2 and nchan > 2:
                            avimages, beamlist, freq_av, crms_av = self.windowaverage_cube(chan_images, rms_desired, chan_rms,
                                        c_wts, sbeam, freqin, nmax_to_avg=nmax_to_avg)
                            total_flux, e_total_flux = self.fit_channels(img, avimages, crms_av, src, beamlist)
                            src_total_flux = N.sum(total_flux, 1) # sum over all Gaussians in source to obtain total fluxes in each channel
                            src_e_total_flux = N.sqrt(N.sum(N.power(e_total_flux, 2.0), 1))
                            src_mask, n_good_chan = self.mask_upper_limits(src_total_flux, src_e_total_flux, snr_desired)
                            nchan = avimages.shape[0]
                            if nchan == nchan_prev:
                                break
                            nchan_prev = nchan
                            rms_desired *= 0.8

                        # Now fit source for spectral index.
                        src_total_flux = src_total_flux.reshape((src_total_flux.shape[0],))
                        src_e_total_flux = src_e_total_flux.reshape((src_e_total_flux.shape[0],))
                        src_mask = src_mask.reshape((src_mask.shape[0],))
                        if img.opts.verbose_fitting:
                            if img.opts.flagchan_snr:
                                print('Source #%i : averaged to %i channels, of which %i meet SNR criterion' % (src.source_id,
                                      len(src_total_flux), nchan))
                            else:
                                print('Source #%i : averaged to %i channels, all of which will be used' % (src.source_id,
                                       len(src_total_flux)))
                        npos = len(N.where(src_total_flux > 0.0)[0])

                        if isinstance(n_good_chan, int):
                            n_good_chan = [n_good_chan]
                        if (img.opts.flagchan_snr and n_good_chan[0] < 2) or npos < 2:
                            src.spec_indx = N.NaN
                            src.e_spec_indx = N.NaN
                            src.spec_norm = N.NaN
                            src.specin_flux = [N.NaN]
                            src.specin_fluxE = [N.NaN]
                            src.specin_freq = [N.NaN]
                            src.specin_freq0 = N.NaN
                        else:
                            if img.opts.flagchan_snr:
                                good_fluxes_ind = N.where(src_mask == False)
                            else:
                                good_fluxes_ind = range(len(freq_av))
                            fluxes_to_fit = src_total_flux[good_fluxes_ind]
                            e_fluxes_to_fit = src_e_total_flux[good_fluxes_ind]
                            freqs_to_fit = freq_av[good_fluxes_ind]

#                             if len(freqs_to_fit.shape) == 2:
#                                 freqs_to_fit = freqs_to_fit.reshape((freqs_to_fit.shape[0],))
#                             if len(fluxes_to_fit.shape) == 2:
#                                 fluxes_to_fit = fluxes_to_fit.reshape((fluxes_to_fit.shape[0],))
#                             if len(e_fluxes_to_fit.shape) == 2:
#                                 e_fluxes_to_fit = e_fluxes_to_fit.reshape((e_fluxes_to_fit.shape[0],))

                            fit_res = self.fit_specindex(freqs_to_fit, fluxes_to_fit, e_fluxes_to_fit)
                            src.spec_norm, src.spec_indx, src.e_spec_indx = fit_res
                            src.specin_flux = fluxes_to_fit.tolist()
                            src.specin_fluxE = e_fluxes_to_fit.tolist()
                            src.specin_freq = freqs_to_fit.tolist()
                            src.specin_freq0 = N.median(freqs_to_fit)
                    else:
                        src.spec_norm = src.gaussians[0].spec_norm
                        src.spec_indx = src.gaussians[0].spec_indx
                        src.e_spec_indx = src.gaussians[0].e_spec_indx
                        src.specin_flux = src.gaussians[0].specin_flux
                        src.specin_fluxE = src.gaussians[0].specin_fluxE
                        src.specin_freq = src.gaussians[0].specin_freq
                        src.specin_freq0 = src.gaussians[0].specin_freq0

                    if bar2.started:
                        bar2.increment()
                if bar2.started:
                    bar2.stop()
                img.completed_Ops.append('spectralindex')
            else:
                mylog.warning('Image has only one channel. Spectral index module disabled.')
                img.opts.spectralindex_do = False

####################################################################################
    def flagchans_rmschan(self, crms, zeroflags, iniflags, cutoff):
        """ Calculate clipped rms (r1) of the rms as fn of channel, crms, with zeroflags
        applied and kappa=cutoff. Then exclude crms=0 (for NaN mages etc) and get ch.s
        which are more than cutoff*r1 away from median of rms. If this is less than 10 %
        of all channels, flag them.

        """
        # crms_rms and median dont include rms=0 channels
        nchan = len(crms)
        mean, rms, cmean, crms_rms, cnt = func.bstat(crms, zeroflags, cutoff)
        zeroind = N.where(crms==0)[0]
        median = N.median(N.delete(crms, zeroind))
        badind = N.where(N.abs(N.delete(crms, zeroind) - median)/crms_rms >=cutoff)[0]
        frac = len(badind)/(nchan - len(zeroind))

        if frac <= 0.1:
            badind = N.where(N.abs(crms - median)/crms_rms >=cutoff)[0]
            iniflags[badind] = True

        return iniflags

####################################################################################
    def iniflag(self, img):
        """ Calculate clipped rms of every channel, and then median and clipped rms of this rms distribution.
        Exclude channels where rms=0 (all pixels 0 or blanked) and of the remaining, if outliers beyond 5 sigma
        are less then 10 % of number of channels, flag them. This is done only when flagchan_rms = True.
        If False, only rms=0 (meaning, entire channel image is zero or blanked) is flagged."""

        image = img.image_arr
        nchan = image.shape[1]
        iniflags = N.zeros(nchan, bool)
        zeroflags = N.zeros(nchan, bool)
        crms = img.channel_clippedrms

        # First, check whether user has specified any channels to flag
        if img.opts.flagchan_list is not None:
            for chan in img.opts.flagchan_list:
                zeroflags[chan] = True

        # Next, flag channels with rms = 0
        for ichan in range(nchan):
            if crms[ichan] == 0: zeroflags[ichan] = True
        iniflags = cp(zeroflags)

        # Lastly, flag outliers
        if img.opts.flagchan_rms:
            iniflags = self.flagchans_rmschan(crms, zeroflags, iniflags, 4.0)

        return iniflags


####################################################################################
    def freq_beamsp_unav(self, img):
        """ Defines img.beam_spectrum and img.freq for the unaveraged cube. """
        # Find the channel frequencies
        shp = img.image_arr.shape
        img.freq = N.zeros(shp[1])
        crval, cdelt, crpix = img.freq_pars
        if img.wcs_obj.wcs.spec == -1 and \
                img.opts.frequency_sp is None:
            raise RuntimeError("Frequency info not found in header "\
                                   "and frequencies not specified by user")
        else:
            if img.opts.frequency_sp is None:
                for ichan in range(shp[1]):
                    img.freq[ichan] = img.wcs_obj.p2f(ichan)
            else:
                if len(img.opts.frequency_sp) != shp[1]:
                    raise RuntimeError("Number of channels does not match number "\
                                 "of frequencies specified by user")
                for ichan in range(shp[1]):
                    img.freq[ichan] = img.opts.frequency_sp[ichan]

        # Find the channel beam shapes
        sbeam = img.opts.beam_spectrum
        if sbeam is not None and len(sbeam) != shp[1]:
            sbeam = None  # sanity check
        if sbeam is None:
            sbeam = []
            hdr = img.header
            try:
                # search for channel beams in the image header
                for ichan in range(shp[1]):
                    sbeam.append((hdr['BMAJ{}'.format(ichan+1)],
                                 hdr['BMIN{}'.format(ichan+1)],
                                 hdr['BPA{}'.format(ichan+1)]))
            except KeyError:
                # Channel beam info not found. Use constant beam or one scaled with
                # frequency
                if img.opts.beam_sp_derive:
                    # Adjust channel beam sizes assuming that the beam scales as 1/nu
                    # Note: beam is (major, minor, pos. angle)
                    for ichan in range(shp[1]):
                        sbeam.append((img.beam[0] * img.freq[0] / img.freq[ichan],
                                     img.beam[1] * img.freq[0] / img.freq[ichan],
                                     img.beam[2]))
                else:
                    sbeam = [img.beam] * shp[1]
        img.beam_spectrum = sbeam

####################################################################################
    def rms_spectrum(self, img, image):
        from .rmsimage import Op_rmsimage
        global bar1
        mylog = img.mylog

        nchan = image.shape[0]
        rms_map = img.use_rms_map
        if img.opts.kappa_clip is None:
            kappa = -img.pixel_beamarea()
        else:
            kappa = img.opts.kappa_clip
        map_opts = (kappa, img.rms_box, img.opts.spline_rank)

        if rms_map:
            rms_spec = N.zeros(image.shape, dtype=N.float32)
            mean = N.zeros(image.shape[1:], dtype=N.float32)
            rms = N.zeros(image.shape[1:], dtype=N.float32)
            median_rms = N.zeros(nchan)
            for ichan in range(nchan):
                if bar1.started:
                    bar1.increment()
                dumi = Op_rmsimage()
                mask = N.isnan(image[ichan])
                Op_rmsimage.map_2d(dumi, image[ichan], mean, rms, mask, *map_opts)
                rms_spec[ichan,:,:] = rms
                median_rms[ichan] = N.median(rms)
        else:
            rms_spec = N.zeros(image.shape, dtype=N.float32)
            for ichan in range(nchan):
                if bar1.started:
                    bar1.increment()
                rms_spec[ichan,:,:] = img.channel_clippedrms[ichan]
            median_rms = rms_spec
        if bar1.started:
            bar1.stop()

        str1 = " ".join(["%9.4e" % n for n in img.channel_clippedrms])
        if rms_map:
            mylog.debug('%s %s ' % ('Median rms of channels : ', str1))
            mylog.info('RMS image made for each channel')
        else:
            mylog.debug('%s %s ' % ('RMS of channels : ', str1))
            mylog.info('Clipped rms calculated for each channel')

        return rms_spec


####################################################################################
    def fit_specindex(self, freqarr, fluxarr, efluxarr, do_log=False):
        """ Fits spectral index to data.

        do_log is True/False implies you fit spectral index in logFlux vs logFreq space or not."""
        from . import functions as func
        import math

        x = freqarr
        flux = fluxarr
        eflux = efluxarr
        f0 = N.median(x)
        mask = N.zeros(len(fluxarr), dtype=bool)
        nan_errors = N.isnan(efluxarr)
        mask[nan_errors] = 1

        if do_log:
            x = N.log10(x/f0); y = N.log10(flux); sig = N.abs(eflux/flux)/2.303
            funct = func.poly
        else:
            x = x/f0; y = flux; sig = eflux
            funct = func.sp_in

        spin, espin = func.fit_mask_1d(x, y, sig, mask, funct, do_err=True, order=1)

        if do_log:
            spin[0] = math.pow(10.0, spin[0])
            espin[0] = spin[0]*math.log(10.0)*espin[0]

        return spin[0], spin[1], espin[1]


########################################################################################

    def windowaverage_cube(self, imagein, rms_desired, chanrms, c_wts, sbeam,
                            freqin, n_min=2, nmax_to_avg=10):
        """Average neighboring channels of cube to obtain desired rms in at least n_min channels

        The clipped rms of each channel is compared to the desired rms. If the
        clipped rms is too high, the channel is averaged with as many neighboring
        channels as necessary to obtain at least the desired rms. This is done
        until the number of OK channels is 2. The averaging is done first at
        the frequency extremes, as the frequency range of the resulting averaged
        flux array will be maximized.

        For example, if the desired rms is 0.1 and the list of rms's is:

            [0.2, 0.2, 0.3, 0.2, 0.2]

        the resulting channels that will be averaged are:

            [[0, 1], [2], [3, 4]]
        """
        from math import sqrt
        from .collapse import avspc_blanks

        # chan_list is a list of lists of channels to average. E.g., if we have
        # 5 channels and we want to average only the first 2:
        #     chan_list = [[0,1], [2], [3], [4]]
        if len(chanrms.shape) ==3:
            crms = N.mean(N.nanmean(chanrms, axis=1), axis=1)
        else:
            crms = chanrms
        chan_list = self.get_avg_chan_list(rms_desired, crms, nmax_to_avg)

        n_new = len(chan_list)
        beamlist = []
        crms_av = N.zeros(n_new)
        freq_av = N.zeros(n_new)
        imageout = N.zeros((n_new, imagein.shape[1], imagein.shape[2]), dtype=N.float32)
        for ichan, avg_list in enumerate(chan_list):
            if len(avg_list) > 1:
                imageout[ichan], dum = avspc_blanks(avg_list, imagein, crms, c_wts)
                chan_slice = slice(avg_list[0], avg_list[1]+1)
                beamlist.append(tuple(N.mean(sbeam[chan_slice], axis=0)))
                freq_av[ichan] = N.mean(freqin[chan_slice])
                crms_av[ichan] = 1.0/sqrt(N.sum(1.0/crms[chan_slice]**2))
            else:
                imageout[ichan] = imagein[avg_list[0]]
                beamlist.append(sbeam[avg_list[0]])
                freq_av[ichan] = N.mean(freqin[avg_list[0]])
                crms_av[ichan] = 1.0/sqrt(N.sum(1.0/crms[avg_list[0]]**2))

        return imageout, beamlist, freq_av, crms_av


    def get_avg_chan_list(self, rms_desired, chanrms, nmax_to_avg):
        """Returns a list of channels to average to obtain given rms_desired
        in at least 2 channels"""
        end = 0
        chan_list = []
        nchan = len(chanrms)
        good_ind = N.where(N.array(chanrms)/rms_desired < 1.0)[0]
        num_good = len(good_ind)
        if num_good < 2:
            # Average channels at start of list
            rms_avg = chanrms[0]
            while rms_avg > rms_desired:
                end += 1
                chan_slice = slice(0, end)
                rms_avg = 1.0/N.sqrt(N.sum(1.0/N.array(chanrms)[chan_slice]**2))
                if end == nchan or end == nmax_to_avg:
                    break
            if end == 0:
                end = 1
            chan_list.append(range(end))
            if end == nchan:
                # This means all channels are averaged into one. If this happens,
                # instead average first half and second half to get two channels
                # and return.
                chan_list = [range(0, int(float(nchan)/2.0)), range(int(float(nchan)/2.0), nchan)]
                return chan_list

            # Average channels at end of list
            rms_avg = chanrms[-1]
            end = nchan
            start = nchan
            while rms_avg > rms_desired:
                start -= 1
                chan_slice = slice(start, end)
                rms_avg = 1.0/N.sqrt(N.sum(1.0/chanrms[chan_slice]/chanrms[chan_slice]))
                if end-start == nmax_to_avg:
                    break

            if start <= max(chan_list[0]):
                # This means we cannot get two averaged channels with desired rms,
                # so just average remaining channels
                chan_list.append(range(max(chan_list[0]), nchan))
            else:
                # First append any channels between those averaged at the start
                # and those at the end
                for i in range(max(chan_list[0])+1, start):
                    chan_list.append([i])
                if start < end:
                    chan_list.append(range(start, end))
        else:
            # No averaging needed
            for i in range(nchan):
                chan_list.append([i])
        return chan_list


    def fit_channels(self, img, chan_images, clip_rms, src, beamlist):
        """Fits normalizations of Gaussians in source to multiple channels

        If unresolved, the size of the Gaussians are adjusted to match the
        channel's beam size (given by beamlist) before fitting.

        Returns array of total fluxes (N_channels x N_Gaussians) and array
        of errors (N_channels x N_Gaussians).
        """
        from . import functions as func
        from .const import fwsig

        isl = img.islands[src.island_id]
        isl_bbox = isl.bbox
        nchan = chan_images.shape[0]
        x, y = N.mgrid[isl_bbox]
        gg = src.gaussians
        fitfix = N.ones(len(gg)) # fit only normalization
        srcmask = isl.mask_active

        total_flux = N.zeros((nchan, len(fitfix))) # array of fluxes: N_channels x N_Gaussians
        errors = N.zeros((nchan, len(fitfix))) # array of fluxes: N_channels x N_Gaussians
        for cind in range(nchan):
            image = chan_images[cind]
            gg_adj = self.adjust_size_by_freq(img.beam, beamlist[cind], gg)
            p, ep = func.fit_mulgaus2d(image, gg_adj, x, y, srcmask, fitfix, adj=True)
            pbeam = img.beam2pix(beamlist[cind])
            bm_pix = (pbeam[0]/fwsig, pbeam[1]/fwsig, pbeam[2])  # IN SIGMA UNITS
            for ig in range(len(fitfix)):
                total_flux[cind, ig] = p[ig*6]*p[ig*6+3]*p[ig*6+4]/(bm_pix[0]*bm_pix[1])
            p = N.insert(p, N.arange(len(fitfix))*6+6, total_flux[cind])
            rms_isl = N.nanmean(clip_rms[cind])
            if N.isnan(rms_isl):
                # If the channel rms is all NaNs, use the average rms value over all
                # channels instead
                rms_isl = N.nanmean(clip_rms)
            if not N.isnan(rms_isl):
                errors[cind] = func.get_errors(img, p, rms_isl, bm_pix=(bm_pix[0]*fwsig, bm_pix[1]*fwsig, bm_pix[2]))[6]
            self.reset_size(gg)

        return total_flux, errors

    def adjust_size_by_freq(self, beam_ch0, beam, gg):
        """Adjust size of unresolved Gaussians to match the channel's beam size"""
        gg_adj = []
        for g in gg:
            g.size_pix_adj = g.size_pix[:]
            if g.deconv_size_sky[0] == 0.0:
                g.size_pix_adj[0] *= beam[0] / beam_ch0[0]
            if g.deconv_size_sky[1] == 0.0:
                g.size_pix_adj[1] *= beam[1] / beam_ch0[1]
            gg_adj.append(g)
        return gg_adj

    def reset_size(self, gg):
        """Reset size of unresolved Gaussians to match the ch0 beam size"""
        for g in gg:
            if hasattr(g, 'size_pix_adj'): del g.size_pix_adj

    def mask_upper_limits(self, total_flux, e_total_flux, threshold):
        """Returns mask of upper limits"""
        mask = N.zeros(total_flux.shape, dtype=bool)
        if len(total_flux.shape) == 1:
            is_src = True
            ndet = 0
            ncomp = 1
        else:
            is_src = False
            ndet = N.zeros((total_flux.shape[1]), dtype=int)
            ncomp = len(ndet)
        for ig in range(ncomp):
            for ichan in range(total_flux.shape[0]):
                if is_src:
                    meas_flux = total_flux[ichan]
                    e_meas_flux = e_total_flux[ichan]
                else:
                    meas_flux = total_flux[ichan, ig]
                    e_meas_flux = e_total_flux[ichan, ig]
                if meas_flux < threshold * e_meas_flux:
                    # Upper limit
                    if is_src:
                        mask[ichan] = True
                    else:
                        mask[ichan, ig] = True
                else:
                    # Detection
                    if is_src:
                        ndet += 1
                        mask[ichan] = False
                    else:
                        ndet[ig] += 1
                        mask[ichan, ig] = False
        return mask, ndet
