.. _process_image:

***********************************************
``process_image``: processing an image
***********************************************

A standard analysis is performed using the ``process_image`` task. This task reads in the input image, calculates background rms and mean images, finds islands of emission, fits Gaussians to the islands, and groups the Gaussians into sources. Furthermore, the ``process_image`` task encompases a number of modules that allow decomposing an image into shapelets, calculating source spectral indices, deriving source polarization properties, and correcting for PSF variations across the image.

When process_image is executed, PyBDSM performs the following steps in
order:

#. Reads in the image and collapses specific frequency channels with weights (see :ref:`multichan_opts`) and produces a 'continuum' image (the ch0 image) for all polarisations with which source detection is done.

#. Calculates basic statistics of the image and sensible values of the processing parameters. First, the number of beams per
   source is calculated (see :ref:`algorithms` for details), using a
   sensible estimate of box size and step size (which can be set using the
   :term:`rms_box` parameter). Next, the thresholds are set. They can either be
   hard thresholded (by the user or set as 5-sigma for pixel threshold and
   3-sigma for island boundaries by default) or can be calculated using the
   False Detection Rate (FDR) method using a user defined value for
   :math:`\alpha` (the :term:`fdr_alpha` parameter). If the user does not specify whether hard thresholding or FDR thresholding
   should be applied, one or the other is chosen internally based on the
   ratio of expected false pixels to true pixels.

#. Calculates the rms and mean images. The 3-sigma clipped rms and mean are calculated
   inside boxes of defined by the :term:`rms_box` parameter. Optionally, these images can be calculated using
   adaptive scaling of this box, so that a smaller box (defined the the :term:`rms_box_bright` parameter) is used near bright sources (where strong artifacts are more likely). Intermediate values
   are calculated using bicubic spline interpolation by default (the order of the spline interpolation can be set with the :term:`spline_rank` parameter). Depending on the resulting statistics (see :ref:`algorithms` for details), we either adopt the rms image or a constant rms
   in the following analysis.

#. Identifies islands of contiguous emission. First all pixels greater
   than the pixel threshold are identified. Next, starting from each of these pixels, all contiguous pixels
   (defined by 8-connectivity, i.e., the surrounding eight pixels) higher
   than the island boundary threshold are identified as belonging to one
   island, accounting properly for overlaps of islands.

#. Fits multiple Gaussians to each island. The number of
   multiple Gaussians to be fit can be determined by three different
   methods (using the :term:`ini_gausfit` parameter). With initial guesses
   corresponding to these peaks, Gaussians are simultaneously fit to the
   island using the Levenberg-Marqhardt algorithm. Sensible criteria for bad
   solutions are defined (see :ref:`flagging_opts`). If multiple Gaussians are fit and one of them is
   a bad solution then the number of Gaussians is decreased by one and fit
   again, until all solutions in the island are good (or zero in number, in
   which case it's flagged). After the final fit to the island, the
   deconvolved size is computed assuming the theoretical beam, and the
   statistics in the source area and in the island are computed and
   stored. Errors on each of the fitted parameters are computed using the
   formulae in Condon (1997) [#f1]_.

#. Groups nearby Gaussians within an island into sources. See :ref:`grouping`
   for details. Fluxes for the grouped Gaussians are summed to obtain the
   total flux of the source (the uncertainty is calculated by summing the
   Gaussian uncertainties in quadrature). The source position is set to be its
   centroid (the position of the maximum of the source is also calculated and
   output). The total source size is measured using moment analysis (see
   http://en.wikipedia.org/wiki/Image_moment for a nice overview of moment
   analysis). Errors on the source position and size are estimated using
   Condon (1997), but additional uncertainties due to uncertainties in the
   constituent Gaussians may be estimated using a Monte Carlo technique.

#. Continues with further processing, if the user has specified that one or more of the following modules be used:

    * Shapelet decomposition (see :ref:`shapelet_do` for details)

    * Wavelet decomposition (see :ref:`atrous_do` for details)

    * Estimation of PSF variation (see :ref:`psf_vary_do` for details)

    * Calculation of polarization properties (see :ref:`polarisation_do` for details)

    * Calculation of spectral indices (see :ref:`spectralindex_do` for details)

.. _general_pars:

General reduction parameters
----------------------------
Type ``inp process_image`` to list the main reduction parameters:

.. parsed-literal::

    PROCESS_IMAGE: Find and measure sources in an image.
    ================================================================================
    :term:`filename` ................. '': Input image file name
    :term:`adaptive_rms_box` ..... False : Use adaptive rms_box when determining rms and
                                   mean maps
    :term:`advanced_opts` ........ False : Show advanced options
    :term:`atrous_do` ............ False : Decompose Gaussian residual image into multiple
                                   scales
    :term:`beam` .................. None : FWHM of restoring beam. Specify as (maj, min, pos
                                   ang E of N) in degrees. E.g., beam = (0.06, 0.02,
                                   13.3). None => get from header
    :term:`flagging_opts` ........ False : Show options for Gaussian flagging
    :term:`frequency` ............. None : Frequency in Hz of input image. E.g., frequency =
                                   74e6. None => get from header. For more than one
                                   channel, use the frequency_sp parameter.
    :term:`interactive` .......... False : Use interactive mode
    :term:`mean_map` .......... 'default': Background mean map: 'default' => calc whether to
                                   use or not, 'zero' => 0, 'const' => clipped mean,
                                   'map' => use 2-D map.
    :term:`multichan_opts` ....... False : Show options for multi-channel images
    :term:`output_opts` .......... False : Show output options
    :term:`polarisation_do` ...... False : Find polarisation properties
    :term:`psf_vary_do` .......... False : Calculate PSF variation across image
    :term:`rms_box` ............... None : Box size, step size for rms/mean map calculation.
                                   Specify as (box, step) in pixels. E.g., rms_box =
                                   (40, 10) => box of 40x40 pixels, step of 10
                                   pixels. None => calculate inside program
    :term:`rms_map` ............... None : Background rms map: True => use 2-D rms map;
                                   False => use constant rms; None => calculate
                                   inside program
    :term:`shapelet_do` .......... False : Decompose islands into shapelets
    :term:`spectralindex_do` ..... False : Calculate spectral indices (for multi-channel
                                   image)
    :term:`thresh` ................ None : Type of thresholding: None => calculate inside
                                   program, 'fdr' => use false detection rate
                                   algorithm, 'hard' => use sigma clipping
    :term:`thresh_isl` ............. 3.0 : Threshold for the island boundary in number of
                                   sigma above the mean. Determines extent of
                                   island used for fitting
    :term:`thresh_pix` ............. 5.0 : Source detection threshold: threshold for the
                                   island peak in number of sigma above the mean. If
                                   false detection rate thresholding is used, this
                                   value is ignored and thresh_pix is calculated
                                   inside the program

Each of the parameters is described in detail below.

.. glossary::
    filename
        This parameter is a string (no default) that sets the input image file name. The input image can be a FITS or CASA 2-, 3-, or 4-D cube.

    adaptive_rms_box
        This parameter is a Boolean (default is ``False``). If ``True``, an adaptive box is used when calculating the rms and mean maps. See :ref:`adaptive_rms_box` for details of the options.

    advanced_opts
        This parameter is a Boolean (default is ``False``). If ``True``, the advanced options are shown. See :ref:`advanced_opts` for details of the advanced options.

    atrous_do
        This parameter is a Boolean (default is ``False``). If ``True``, wavelet decomposition will be performed. See :ref:`atrous_do` for details of the options.

    beam
        This parameter is a tuple (default is ``None``) that defines the FWHM of restoring beam. Specify as (maj, min, pos ang E of N) in degrees. E.g., ``beam = (0.06, 0.02, 13.3)``. For more than one channel, use the ``beam_spectrum`` parameter. If the beam is not given by the user, then it is looked for in the image header. If not found, then an error is raised. PyBDSM will not work without knowledge of the restoring beam.

    flagging_opts
        This parameter is a Boolean (default is ``False``). If ``True``, the Gaussian flagging options will be listed. See :ref:`flagging_opts` for details of the options.

    frequency
        This parameter is a float (default is ``None``) that defines the frequency in Hz of the input image. E.g., ``frequency = 74e6``. For more than one channel, use the :term:`frequency_sp` parameter. If the frequency is not given by the user, then it is looked for in the image header. If not found, then an error is raised. PyBDSM will not work without knowledge of the frequency.

    interactive
        This parameter is a Boolean (default is ``False``). If ``True``, interactive mode is used. In interactive mode, plots are displayed at various stages of the processing so that the user may check the progress of the fit.

        First, plots of the rms and mean background images are displayed along with the islands found, before fitting of Gaussians takes place. The user should verify that the islands and maps are reasonable before preceding.

        Next, if ``atrous_do = True``, the fits to each wavelet scale are shown. The wavelet fitting may be truncated at the current scale if desired.

        Lastly, the final results are shown.

    mean_map
        This parameter is a string (default is ``'default'``) that determines how the background mean map is computed and
        how it is used further.

        If ``'const'``\, then the value of the clipped mean of the entire image (set
        by the ``kappa_clip`` option) is used as the background mean map.

        If ``'zero'``\, then a value of zero is used.

        If ``'map'``\, then the 2-dimensional mean map is computed and used. The
        resulting mean map is largely determined by the value of the ``rms_box``
        parameter (see the ``rms_box`` parameter for more information).

        If ``'default'``\, then PyBDSM will attempt to determine automatically
        whether to use a 2-dimensional map or a constant one as follows. First,
        the image is assumed to be confused if ``bmpersrc_th`` < 25 or the ratio of
        the clipped mean to rms (clipped mean/clipped rms) is > 0.1, else the
        image is not confused. Next, the mean map is checked to see if its
        spatial variation is significant. If so, then a 2-D map is used and, if
        not, then the mean map is set to either 0.0 or a constant depending on
        whether the image is thought to be confused or not.

        Generally, ``'default'`` works well. However, if there is significant
        extended emission in the image, it is often necessary to force the use
        of a constant mean map using either ``'const'`` or ``'mean'``\.

    multichan_opts
        This parameter is a Boolean (default is ``False``). If ``True``, the multichannel options will be listed. See :ref:`multichan_opts` for details of the options.

    output_opts
        This parameter is a Boolean (default is ``False``). If ``True``, the output options will be listed. See :ref:`output_opts` for details of the options.

    polarisation_do
        This parameter is a Boolean (default is ``False``). If ``True``, polarization properties will be calculated for the sources. See :ref:`polarisation_do` for details of the options.

    psf_vary_do
        This parameter is a Boolean (default is ``False``). If ``True``, the spatial variation of the PSF will be estimated and its effects corrected. See :ref:`psf_vary_do` for details of the options.

    rms_box
        This parameter is a tuple (default is ``None``) of two integers and is probably the most important input
        parameter for PyBDSM. The first integer, boxsize, is the size of the 2-D
        sliding box for calculating the rms and mean over the entire image. The
        second, stepsize, is the number of pixels by which this box is moved for
        the next measurement. If ``None``\, then suitable values are calculated
        internally.

        In general, it is best to choose a box size that corresponds to the
        typical scale of artifacts in the image, such as those that are common
        around bright sources. Too small of a box size will effectively raise
        the local rms near a source so much that a source may not be fit at all;
        too large a box size can result in underestimates of the rms due to
        oversmoothing. A step size of 1/3 to 1/4 of the box size usually works
        well.

        .. note::

            The :term:`spline_rank` parameter also affects the rms and mean maps. If you find ringing artifacts in the rms or mean maps near bright sources, try adjusting this parameter.

    rms_map
        This parameter is a Boolean (default is ``None``). If ``True``\, then the 2-D background rms image is computed and used. If
        ``False``\, then a constant value is assumed (use ``rms_value`` to force the rms
        to a specific value). If ``None``\, then the 2-D rms image is calculated, and
        if the variation is statistically significant then it is taken, else a
        constant value is assumed. The rms image used for each channel in
        computing the spectral index follows what was done for the
        channel-collapsed image.

        Generally, the default value works well. However, if there is significant extended
        emission in the image, it is often necessary to force the use of a
        constant rms map by setting ``rms_map = False``.

    shapelet_do
        This parameter is a Boolean (default is ``False``). If ``True``, shapelet decomposition of the islands will be performed. See :ref:`shapelet_do` for details of the options.

    spectralindex_do
        This parameter is a Boolean (default is ``False``). If ``True``, spectral indices will be calculated for the sources. See :ref:`spectralindex_do` for details of the options.

    thresh
        This parameter is a string (default is ``None``). If ``thresh = 'hard'``\, then a hard threshold is assumed, given by
        thresh_pix. If ``thresh = 'fdr'``\, then the False Detection Rate algorithm
        of Hopkins et al. (2002) is used to calculate the value of ``thresh_pix``\.
        If ``thresh = None``\, then the false detection probability is first
        calculated, and if the number of false source pixels is more than
        ``fdr_ratio`` times the estimated number of true source pixels, then the
        ``'fdr'`` threshold option is chosen, else the ``'hard'`` threshold option is
        chosen.

    thresh_isl
        This parameter is a float (default is 3.0) that determines the region to which fitting is done. A higher
        value will produce smaller islands, and hence smaller regions that are
        considered in the fits. A lower value will produce larger islands. Use
        the thresh_pix parameter to set the detection threshold for sources.
        Generally, ``thresh_isl`` should be lower than ``thresh_pix``\.

        Only regions above the absolute threshold will be used. The absolute
        threshold is calculated as ``abs_thr = mean + thresh_isl * rms``\. Use the
        ``mean_map`` and ``rms_map`` parameters to control the way the mean and rms are
        determined.

    thresh_pix
        This parameter is a float (default is 5.0) that sets the source detection threshold in number of
        sigma above the mean. If false detection rate thresholding is used, this
        value is ignored and ``thresh_pix`` is calculated inside the program

        This parameter sets the overall detection threshold for islands (i.e.
        ``thresh_pix = 5`` will find all sources with peak flux densities per beam of 5-sigma or
        greater). Use the ``thresh_isl`` parameter to control how much of each
        island is used in fitting. Generally, ``thresh_pix`` should be larger than
        ``thresh_isl``.

        Only islands with peaks above the absolute threshold will be used. The
        absolute threshold is calculated as ``abs_thr = mean + thresh_pix * rms``\.
        Use the ``mean_map`` and ``rms_map`` parameters to control the way the mean and
        rms are determined.


.. _adaptive_rms_box:

Adaptive box options
====================
If ``adaptive_rms_box = True``, the rms_box is reduced in size near bright sources and enlarged far from them. This scaling attempts to account for possible strong artifacts around bright sources while still acheiving accurate background rms and mean values when extended sources are present. This option is generally slower than non-adaptive scaling.

Use the ``rms_box`` parameter to set the large-scale box and the ``rms_box_bright`` parameter to set the small-scale box. The threshold for bright sources can be set with the ``adaptive_thresh`` parameter:

.. parsed-literal::

    adaptive_rms_box ...... True : Use adaptive rms_box when determining rms and mean maps
      :term:`adaptive_thresh` ..... None : Sources with pixels above adaptive_thresh*
                                   clipped_rms will be considered as bright sources (i.e.,
                                   with potential artifacts). None => calculate inside
                                   program
      :term:`rms_box_bright` ...... None : Box size, step size for rms/mean map
                                   calculation near bright sources. Specify as (box, step)
                                   in pixels. None => calculate inside program

.. glossary::

    adaptive_thresh
        This parameter is a float (default is ``None``) that sets the SNR above which sources may be affected by strong artifacts Sources that meet the SNR threshold will use the small-scale box (set by the ``rms_box_bright`` parameter) if their sizes at a threshold of 10.0 is less than 25 beam areas.

        If None, the threshold is varied from 500 to 50 to attempt to obtain at least 5 candidate bright sources.

    rms_box_bright
        This parameter is a tuple (default is ``None``) of two integers that sets the box and step sizes to use near bright sources (determined by the ``adaptive_thresh`` parameter). The large-scale box size is set with the ``rms_box`` parameter.

.. _advanced_opts:

Advanced options
================
If ``advanced_opts = True``, a number of additional options are listed. The advanced options do not usually need to be altered from the default values, but can be useful, for example, for fine tuning a fit or for quickly fitting a small region of a much larger image.

The advanced options are:

.. parsed-literal::

    advanced_opts ......... True : Show advanced options
      :term:`aperture` ............ None : Radius of aperture in pixels inside which aperture
                                   fluxes are measured for each source. None => no aperture
                                   fluxes measured
      :term:`blank_zeros` ........ False : Blank zeros in the image
      :term:`bmpersrc_th` ......... None : Theoretical estimate of number of beams per
                                   source. None => calculate inside program
      :term:`check_outsideuniv` .. False : Check for pixels outside the universe
      :term:`detection_image` ........ '': Detection image file name used only for
                                   detecting islands of emission. Source
                                   measurement is still done on the main image
      :term:`do_mc_errors` ....... False : Estimate uncertainties for 'M'-type sources
                                   using Monte Carlo method
      :term:`fdr_alpha` ........... 0.05 : Alpha for FDR algorithm for thresholds
      :term:`fdr_ratio` ............ 0.1 : For thresh = None; if #false_pix / #source_pix <
                                   fdr_ratio, thresh = 'hard' else thresh = 'fdr'
      :term:`fittedimage_clip` ..... 0.1 : Sigma for clipping Gaussians while creating fitted
                                   image
      :term:`group_by_isl` ....... False : Group all Gaussians in each island into a single
                                   source
      :term:`group_tol` ............ 1.0 : Tolerance for grouping of Gaussians into sources:
                                   larger values will result in larger sources
      :term:`ini_gausfit` ..... 'default': Initial guess for Gaussian parameters: 'default',
                                   'fbdsm', or 'nobeam'
      :term:`kappa_clip` ........... 3.0 : Kappa for clipped mean and rms
      :term:`minpix_isl` .......... None : Minimal number of pixels with emission per island.
                                   None -> calculate inside program
      :term:`ncores` .............. None : Number of cores to use during fitting, None => use
                                   all
      :term:`peak_fit` ............ True : Find and fit peaks of large islands before fitting
                                   entire island
      :term:`peak_maxsize` ........ 30.0 : If island size in beam area is more than this,
                                   attempt to fit peaks separately (if
                                   peak_fit=True). Min value is 30
      :term:`rms_value` ........... None : Value of constant rms in Jy/beam to use if rms_map
                                   = False. None => calculate inside program
      :term:`spline_rank` ............ 3 : Rank of the interpolating function for rms/mean
                                   map
      :term:`split_isl` ........... True : Split island if it is too large, has a large
                                   convex deficiency and it opens well. If it doesn't
                                   open well, then isl.mean = isl.clipped_mean, and
                                   is taken for fitting. Splitting, if needed, is
                                   always done for wavelet images
      :term:`splitisl_maxsize` .... 50.0 : If island size in beam area is more than this,
                                   consider splitting island. Min value is 50
      :term:`stop_at` ............. None : Stops after: 'isl' = island finding step or 'read'
                                   = image reading step
      :term:`trim_box` ............ None : Do source detection on only a part of the image.
                                   Specify as (xmin, xmax, ymin, ymax) in pixels.
                                   E.g., trim_box = (120, 840, 15, 895). None => use
                                   entire image

.. glossary::

    aperture
        This parameter is a float (default is ``None``) that sets the radius (in
        pixels) inside which the aperture flux is measured for each source.
        The aperture is centered on the centroid of the source. Errors are
        calculated from the mean of the rms map inside the aperture.

    blank_zeros
        This parameter is a Boolean (default is ``False``). If ``True``, all
        pixels in the input image with values of 0.0 are blanked. If ``False``,
        any such pixels are left unblanked (and hence will affect the rms and
        mean maps, etc.). Pixels with a value of NaN are always blanked.

    bmpersrc_th
        This parameter is a float (default is ``None``) that sets the
        theoretical estimate of number of beams per source. If ``None``, the
        value is calculated as N/[n*(alpha-1)], where N is the total number of
        pixels in the image, n is the number of pixels in the image whose value
        is greater than 5 times the clipped rms, and alpha is the slope of the
        differential source counts distribution, assumed to be 2.5.

        The value of ``bmpersrc_th`` is used
        to estimate the average separation in pixels between two sources, which
        in turn is used to estimate the boxsize for calculating the background
        rms and mean images. In addition, if the value is below 25 (or the ratio
        of clipped mean to clipped rms of the image is greater than 0.1), the
        image is assumed to be confused and hence the background mean is put to
        zero.

    check_outsideuniv
        This parameter is a Boolean (default is ``False``). If ``True``, then
        the coordinate of each pixel is examined to check if it is outside the
        universe, which may happen when, e.g., an all sky image is made with SIN
        projection (commonly done at LOFAR earlier). When found, these pixels
        are blanked (since imaging software do not do this on their own). Note
        that this process takes a lot of time, as every pixel is checked in case
        weird geometries and projections are used.

    detection_image
        This parameter is a string (default is ``''``) that sets the detection
        image file name used only for detecting islands of emission. Source
        measurement is still done on the main image. The detection image can be
        a FITS or CASA 2-, 3-, or 4-D cube and must have the same size and WCS
        parameters as the main image.

    do_mc_errors
        This parameter is a Boolean (default is ``False``). If ``True``,
        uncertainties on the sizes and positions of 'M'-type sources due to
        uncertainties in the constituent Gaussians are estimated using a Monte
        Carlo technique. These uncertainties are added in quadrature with those
        calculated using Condon (1997). If ``False``, these uncertainties are
        ignored, and errors are calculated using Condon (1997) only.

        Enabling this option will result in longer run times if many 'M'-type
        sources are present, but should give better estimates of the
        uncertainites, particularly for complex sources composed of many
        Gaussians.

    fdr_alpha
        This parameter is a float (default is 0.05) that sets the value of alpha
        for the FDR algorithm for thresholding. If ``thresh`` is ``'fdr'``, then
        the estimate of ``fdr_alpha`` (see Hopkins et al. 2002 [#f2]_ for
        details) is stored in this parameter.

    fdr_ratio
        This parameter is a float (default is 0.1). When ``thresh = None``, if
        #false_pix / #source_pix < fdr_ratio, ``thresh = 'hard'`` otherwise
        ``thresh = 'fdr'``.

    fittedimage_clip
        This parameter is a float (default is 0.1). When the residual image is
        being made after Gaussian decomposition, the model images for each
        fitted Gaussian are constructed up to a size 2b, such that the amplitude
        of the Gaussian falls to a value of ``fitted_image_clip`` times the
        local rms, b pixels from the peak.

    group_by_isl
        This parameter is a Boolean (default is ``False``). If True, all
        Gaussians in the island belong to a single source. If False, grouping is
        controlled by the group_tol parameter.

    group_tol
        This parameter is a float (default is 1.0) that sets the tolerance for grouping of Gaussians into sources: larger values will
        result in larger sources. Sources are created by grouping nearby Gaussians as follows: (1) If the
        minimum value between two Gaussians in an island is more than ``group_tol * thresh_isl * rms_clip``\, and (2) if the centres are seperated by a
        distance less than 0.5*``group_tol`` of the sum of their FWHMs along the PA
        of the line joining them, they belong to the same island.

    ini_gausfit
        This parameter is a string (default is ``'default'``). These are three different ways of estimating the initial guess for
        fitting of Gaussians to an island of emission. If ``'default'``, the maximum
        number of Gaussians is estimated from the number of peaks in the island.
        An initial guess is made for the parameters of these Gaussians before
        final fitting is done. This method should produce the best results when
        there are no large sources present. If ``'simple'``, the maximum number of
        Gaussians per island is set to 25, and no initial guess for the Gaussian
        parameters is made. Lastly, the ``'nobeam'`` method is similar to the
        ``'default'`` method, but no information about the beam is used. This method
        is best used when source sizes are expected to be very different from
        the beam and is generally slower than the other methods. For wavelet
        images, the value used for the original image is used for wavelet order
        j <= 3 and ``'nobeam'`` for higher orders.

    kappa_clip
        This parameter is a float (default is 3.0) that is the factor used for Kappa-alpha clipping, as in
        AIPS. For an image with few source pixels added on to (Gaussian) noise
        pixels, the dispersion of the underlying noise will need to be
        determined. This is done iteratively, whereby the actual dispersion is
        first computed. Then, all pixels whose value exceeds kappa clip times
        this rms are excluded and the rms is computed again. This process is
        repeated until no more pixels are excluded. For well behaved noise
        statistics, this process will converge to the true noise rms with a
        value for this parameter ~3-5. A large fraction of source pixels, less
        number of pixels in total, or significant non-Gaussianity of the
        underlying noise will all lead to non-convergence.

    minpix_isl
        This parameter is an integer (default is ``None``) that sets the minimum number of pixels in an island
        for the island to be included. If
        ``None``\, the number of pixels is set to 1/3 of the area of an unresolved source
        using the beam and pixel size information in the image header. It is set
        to 6 pixels for all wavelet images.

    ncores
        This parameter is an integer (default is ``None``) that sets the number of cores to use during fitting.

    peak_fit
        This parameter is a Boolean (default is ``True``). When True, PyBDSM will identify and fit peaks of emission in large islands iteratively (the size of islands for which peak fitting is done is controlled with the peak_maxsize option), using a maximum of 10 Gaussians per iteration. Enabling this option will generally speed up fitting (by factors of many for large islands), but may result in somewhat higher residuals.

    peak_maxsize
        This parameter is a float (default is 30.0). If island size in beam area is more than this value, attempt to fit peaks
        iteratively (if ``peak_fit = True``). The minimum value is 30.

    rms_value
        This parameter is a float (default is ``None``) that sets the value of constant rms in Jy/beam to use if ``rms_map = False``. If ``None``, the value is
        calculated inside the program.

    spline_rank
        This parameter is an integer (default is 3) that sets the order of the interpolating spline function
        to interpolate the background rms and mean maps over the entire image.

        .. note::

            Bicubic interpolation (the default) can cause ringing artifacts to appear in the rms and mean maps in regions where sharp changes occur. If you find such artifacts, try changing the :term:`spline_rank` parameter.

    split_isl
        This parameter is a Boolean (default is ``True``). If ``True``, an island is split if it is too large, has a large convex deficiency and it
        opens well. If it doesn't open well, then ``isl.mean = isl.clipped_mean``,
        and is taken for fitting. Splitting, if needed, is always done for
        wavelet images

    splitisl_maxsize
        This parameter is a float (default is 50.0). If island size in beam area is more than this, consider splitting
        island. The minimum value is 50.

    stop_at
        This parameter is a string (default is ``None``) that stops an analysis after: 'isl' = island finding step or 'read' = image reading step.

    trim_box
        This parameter is a tuple (default is ``None``) that defines a subregion of the image on which to do source detection. It is specified as (xmin, xmax,
        ymin, ymax) in pixels. E.g., ``trim_box = (120, 840, 15, 895)``\. If ``None``, the entire image is used.


.. _flagging_opts:

Flagging options
================
If ``flagging_opts = True``, a number of options are listed for flagging unwanted Gaussians that occur durring a fit. Flagged Gaussians are not included in any further analysis or catalog. They may be viewed using the ``show_fit`` task (see :ref:`showfit`). A flag value is associated with each flagged Gaussian that allows the user to determine the reason or reasons that it was flagged. If multiple flagging conditions are met by a single Gaussian, the flag values are summed. For example, if a Gaussian is flagged because it is too large (its size exceeds that implied by ``flag_maxsize_bm``, giving a flag value of 64) and because it is too bright (its peak flux density per beam exceeds that implied by ``flag_maxsnr``, giving a flag value of 2) then the final flag value is 64 + 2 = 66.

.. note::

    If a fit did not produce good results, it is often useful to check whether there are flagged Gaussians and adjust the flagging options as necessary. Flagged Gaussians can be viewed by setting ``ch0_flagged = True`` in the ``show_fit`` task.

The options for flagging of Gaussians are:

.. parsed-literal::

    flagging_opts ......... True : Show options for Gaussian flagging
      :term:`flag_bordersize` ........ 0 : Flag Gaussian if centre is outside border -
                                   flag_bordersize pixels
      :term:`flag_maxsize_bm` ..... 50.0 : Flag Gaussian if area greater than flag_maxsize_bm
                                   times beam area
      :term:`flag_maxsize_isl` ..... 1.0 : Flag Gaussian if x, y bounding box around
                                   sigma-contour is factor times island bbox
      :term:`flag_maxsnr` .......... 1.5 : Flag Gaussian if peak is greater than flag_maxsnr
                                   times max value in island
      :term:`flag_minsize_bm` ...... 0.7 : Flag Gaussian if flag_smallsrc = True and area
                                   smaller than flag_minsize_bm times beam area
      :term:`flag_minsnr` .......... 0.9 : Flag Gaussian if peak is less than flag_minsnr
                                   times thresh_pix times local rms
      :term:`flag_smallsrc` ...... False : Flag sources smaller than flag_minsize_bm times
                                   beam area

.. glossary::

    flag_bordersize
        This parameter is an integer (default is 0). Any fitted Gaussian whose centre is ``flag_bordersize`` pixels outside the island
        bounding box is flagged. The flag value is increased by 4 (for x) and 8
        (for y).

    flag_maxsize_bm
        This parameter is a float (default is 25.0). Any fitted Gaussian whose size is greater than ``flag_maxsize_bm`` times the
        synthesized beam is flagged. The flag value is increased by 64.

    flag_maxsize_fwhm
        This parameter is a float (default is 0.3). Any fitted Gaussian whose contour of ``flag_maxsize_fwhm`` times the FWHM falls outside the island is flagged. The flag value is increased by 256.

    flag_maxsize_isl
        This parameter is a float (default is 1.0). Any fitted Gaussian whose maximum x-dimension is larger than
        ``flag_maxsize_isl`` times the x-dimension of the island (and likewise for
        the y-dimension) is flagged. The flag value is increased by 16 (for x)
        and 32 (for y).

    flag_maxsnr
        This parameter is a float (default is 1.5). Any fitted Gaussian whose peak is greater than ``flag_maxsnr`` times
        the value of the image at the peak of the Gaussian is flagged. The flag value is increased
        by 2.

    flag_minsize_bm
        This parameter is a float (default is 0.7). If ``flag_smallsrc`` is True, then any fitted Gaussian whose size is less
        than ``flag_maxsize_bm`` times the synthesized beam is flagged. The Gaussian
        flag is increased by 128.

    flag_minsnr
        This parameter is a float (default is 0.7). Any fitted Gaussian whose peak is less than ``flag_minsnr`` times ``thresh_pix``
        times the local rms is flagged. The flag value is increased by 1.

    flag_smallsrc
        This parameter is a Boolean (default is ``False``). If ``True``\, then fitted Gaussians whose size is less than ``flag_minsize_bm``
        times the synthesized beam area are flagged.  When combining Gaussians
        into sources, an error is raised if a 2x2 box with the peak of the
        Gaussian does not have all four pixels belonging to the source. Usually
        this means that the Gaussian is an artifact or has a very small size.

        If ``False``\, then if either of the sizes of the fitted Gaussian is zero,
        then the Gaussian is flagged.

        If the image is barely Nyquist sampled, this flag is best set to ``False``\.
        This flag is automatically set to ``False`` while decomposing wavelet images
        into Gaussians.

.. _output_opts:

Output options
==============
If ``output_opts = True``, options to control the output generated by ``process_image`` are listed. By default, only a log file is generated and output is controlled with the ``export_image`` (see :ref:`export_image`) and ``write_catalog`` (see :ref:`write_catalog`) tasks. However, the user can specify that a number of optional output files be made automatically whenever ``process_image`` is run. These options are most useful for debugging or when running PyBDSM non-interactively in a Python script (see :ref:`scripting`).

The output options are:

.. parsed-literal::

    output_opts ........... True : Show output options
      :term:`bbs_patches` ......... None : For BBS format, type of patch to use: None => no
                                   patches. 'single' => all Gaussians in one patch.
                                   'gaussian' => each Gaussian gets its own patch.
                                   'source' => all Gaussians belonging to a single
                                   source are grouped into one patch
      :term:`indir` ............... None : Directory of input FITS files. None => get from
                                   filename
      :term:`opdir_overwrite` .. 'overwrite': 'overwrite'/'append': If output_all=True,
                                   delete existing files or append a new directory
      :term:`output_all` ......... False : Write out all files automatically to directory
                                   'filename_pybdsm'
      :term:`output_fbdsm` ....... False : write out fBDSM format output files for use in
                                   Anaamika
      :term:`plot_allgaus` ....... False : Make a plot of all Gaussians at the end
      :term:`plot_islands` ....... False : Make separate plots of each island during fitting
                                   (for large images, this may take a long time and a
                                   lot of memory)
      :term:`print_timing` ....... False : Print basic timing information
      :term:`quiet` .............. False : Suppress text output to screen. Output is still
                                   sent to the log file as usual
      :term:`savefits_meanim` .... False : Save background mean image as fits file
      :term:`savefits_normim` .... False : Save norm image as fits file
      :term:`savefits_rankim` .... False : Save island rank image as fits file
      :term:`savefits_residim` ... False : Save residual image as fits file
      :term:`savefits_rmsim` ..... False : Save background rms image as fits file
      :term:`solnname` ............ None : Name of the run, to be appended to the name of the
                                   output directory
      :term:`verbose_fitting` .... False : Print out extra information during fitting

.. glossary::

    bbs_patches
        This parameter is a string (default is ``None``) that sets the type of patch to use in BBS-formatted catalogs. When the Gaussian catalogue is written as a BBS-readable sky file, this
        determines whether all Gaussians are in a single patch (``'single'``), there are no
        patches (``None``), all Gaussians for a given source are in a separate patch (``'source'``), or
        each Gaussian gets its own patch (``'gaussian'``).

        If you wish to have patches defined by island, then set
        ``group_by_isl = True`` before fitting to force all
        Gaussians in an island to be in a single source. Then set
        ``bbs_patches = 'source'`` when writing the catalog.

    indir
        This parameter is a string (default is ``None``) that sets the directory of input FITS files. If ``None``, the directory is defined by the input filename.

    opdir_overwrite
        This parameter is a string (default is ``'overwrite'``) that determines whether existing output files are overwritten or not.

    output_all
        This parameter is a Boolean (default is ``False``). If ``True``\, all output products are written automatically to the directory ``'filename_pybdsm'``.

    output_fbdsm
        This parameter is a Boolean (default is ``False``). If ``True``\, write out fBDSM format output files for use in Anaamika.

    plot_allgaus
        This parameter is a Boolean (default is ``False``). If ``True``\, make a plot of all Gaussians at the end.

    plot_islands
        This parameter is a Boolean (default is ``False``). If ``True``\, make separate plots of each island during fitting
        (for large images, this may take a long time and a
        lot of memory).

    print_timing
        This parameter is a Boolean (default is ``False``). If ``True``\, print basic timing information.

    quiet
        This parameter is a Boolean (default is ``False``). If ``True``\, suppress text output to screen. Output is still
        sent to the log file as usual.

    savefits_meanim
        This parameter is a Boolean (default is ``False``). If ``True``\, save background mean image as a FITS file.

    savefits_normim
        This parameter is a Boolean (default is ``False``). If ``True``\, save norm image as a FITS file.

    savefits_rankim
        This parameter is a Boolean (default is ``False``). If ``True``\, save island rank image as a FITS file.

    savefits_residim
        This parameter is a Boolean (default is ``False``). If ``True``\, save residual image as a FITS file.

    savefits_rmsim
        This parameter is a Boolean (default is ``False``). If ``True``\, save background rms image as a FITS file.

    solnname
        This parameter is a string (default is ``None``) that sets the name of the run, to be appended to the name of the
        output directory.

    verbose_fitting
        This parameter is a Boolean (default is ``False``). If ``True``\, print out extra information during fitting.



.. _multichan_opts:

Multichannel options
====================
If ``multichan_opts = True``, the options used to control the way PyBDSM handles images with more than one frequency channel are listed. In particular, these options control how the multichannel image is collapsed to form the ``ch0`` image on which source detection is done.

The options concerning multichannel images are:

.. parsed-literal::

    multichan_opts ........ True : Show options for multi-channel images
      :term:`beam_sp_derive` ..... False : If True and beam_spectrum is None, then assume
                                   header beam is for median frequency and scales
                                   with frequency for channels
      :term:`beam_spectrum` ....... None : FWHM of synthesized beam per channel. Specify as
                                   [(bmaj_ch1, bmin_ch1, bpa_ch1), (bmaj_ch2,
                                   bmin_ch2, bpa_ch2), etc.] in degrees. E.g.,
                                   beam_spectrum = [(0.01, 0.01, 45.0), (0.02, 0.01,
                                   34.0)] for two channels. None => all equal to beam
      :term:`collapse_av` ........... [] : List of channels to average if collapse_mode =
                                   'average'; None => all
      :term:`collapse_ch0` ........... 0 : Number of the channel for source extraction, if
                                   collapse_mode = 'single'
      :term:`collapse_mode` ... 'average': Collapse method: 'average' or 'single'. Average
                                   channels or take single channel to perform source
                                   detection on
      :term:`collapse_wt` ....... 'unity': Weighting: 'unity' or 'rms'. Average channels with
                                   weights = 1 or 1/rms_clip^2 if collapse_mode =
                                   'average'
      :term:`frequency_sp` ........ None : Frequency in Hz of channels in input image when
                                   more than one channel is present. E.g., frequency
                                   = [74e6, 153e6]. None => get from header

.. glossary::

    beam_sp_derive
        This parameter is a Boolean (default is ``False``). If `True` and the parameter beam_spectrum is ``None``, then we assume that the
        beam in the header is for the median frequency of the image cube and
        scale accordingly to calculate the beam per channel. If ``False``, then a
        constant value of the beam is taken instead.

    beam_spectrum
        his parameter is a list of tuples (default is ``None``) that sets the FWHM of synthesized beam per channel. Specify as [(bmaj_ch1, bmin_ch1,
        bpa_ch1), (bmaj_ch2, bmin_ch2, bpa_ch2), etc.] in degrees. E.g.,
        ``beam_spectrum = [(0.01, 0.01, 45.0), (0.02, 0.01, 34.0)]`` for two
        channels.

        If ``None``, then the channel-dependent restoring beam is either assumed to
        be a constant or to scale with frequency, depending on whether the
        parameter ``beam_sp_derive`` is ``False`` or ``True``.

    collapse_av
        This parameter is a list of integers (default is ``[]``) that specifies the channels to be averaged to produce the
        continuum image for performing source extraction, if ``collapse_mode`` is
        ``'average'``. If the value is ``[]``, then all channels are used. Otherwise, the
        value is a Python list of channel numbers.

    collapse_ch0
        This parameter is an integer (default is 0) that specifies the number of the channel for source extraction, if ``collapse_mode = 'single'``.

    collapse_mode
        This parameter is a string (default is ``'average'``) that determines whether, when multiple channels are present,
        the source extraction is done on a single channel (``'single'``) or an average of many
        channels (``'average'``).

    collapse_wt
        This parameter is a string (default is ``'unity'``). When ``collapse_mode`` is ``'average'``, then if this value is ``'unity'``, the
        channels given by ``collapse_av`` are averaged with unit weights and if
        ``'rms'``, then they are averaged with weights which are inverse square of
        the clipped rms of each channel image.

    frequency_sp
        This parameter is a list of floats (default is ``None``) that sets the frequency in Hz of channels in input image when more than one channel is present. E.g., ``frequency_sp = [74e6, 153e6]``.

        If the frequency is not given by the user, then it is looked for in the
        image header. If not found, then an error is raised. PyBDSM will not
        work without the knowledge of the frequency.


.. _atrous_do:

*À trous* wavelet decomposition module
--------------------------------------
If ``atrous_do = True``, this module decomposes the residual image that results from the normal fitting of Gaussians into wavelet images of various scales. Such a decomposition is useful if there is extended emission that is not well fit during normal fitting. Such emission therefore remains in the Gaussian residual image and can be further fit by Gaussians whose size is tuned to the various wavelet scales. Therefore, wavelet decomposition should be used when there is significant residual emission that remains after normal Gaussian fitting.

The wavelet module performs the following steps:

* The number of wavelet scales to be considered is set by the ``atrous_jmax`` parameter. By default, this number is determined automatically from the size of the largest island in the image. Wavelet images are then made for scales of order (*j*) ranging from 1 to *jmax*.

* For each scale (*j*), the appropriate *à trous* wavelet transformation is made (see Holschneider et al. 1989 for details). Additionally, the "remainder" image (called the *c_J* image) is also made. This image includes all emission not included in the other wavelet images.

* If ``atrous_bdsm = True``, an rms map and a mean map are made for each wavelet image and Gaussians are fit in the normal way. These wavelet Gaussians can then be included in source catalogs (see :ref:`write_catalog`).

The options for this module are as follows:

.. parsed-literal::

    atrous_do ............. True : Decompose Gaussian residual image into multiple
                                   scales
      :term:`atrous_bdsm_do` ...... True : Perform source extraction on each wavelet scale
      :term:`atrous_jmax` ............ 0 : Max allowed wavelength order, 0 => calculate
                                   inside program
      :term:`atrous_lpf` ........... 'b3': Low pass filter, either 'b3' or 'tr', for B3
                                   spline or Triangle

.. glossary::

    atrous_bdsm_do
        This parameter is a Boolean (default is ``False``). If ``True``\, PyBDSM performs source extraction on each wavelet scale.

    atrous_jmax
        This parameter is an integer (default is 0) which is the maximum order of the *à trous* wavelet
        decomposition. If 0 (or <0 or >15), then the value is determined within
        the program. The value of this parameter is then estimated as the
        (lower) rounded off value of ln[(nm-l)/(l-1) + 1]/ln2 + 1 where nm is
        the minimum of the residual image size (n, m) in pixels and l is the
        length of the filter *à trous* lpf (see the ``atrous_lpf`` parameter for more
        info).

        A sensible value is such that the size of the kernel is not more than
        3-4 times smaller than the smallest image dimension.

    atrous_lpf
        This parameter is a string (default is ``'b3'``) that sets the low pass filter, which can currently be either the B3 spline
        or the triangle function, which is used to generate the *à trous*
        wavelets. The B3 spline is [1, 4, 6, 4, 1] and the triangle is [1, 2,
        1], normalised so that the sum is unity. The lengths of the filters are
        hence 5 and 3 respectively.

.. _psf_vary_do:

PSF variation module
--------------------
If ``psf_vary_do = True``, then the spatial variations in the PSF are estimated and their effects corrected for. To this end, PyBDSM performs the following steps:

* A list of sources that are likely to be unresolved is constructed. This is done by first selecting only type 'S' sources by default (see :ref:`output_cols` for details of source types), but this restriction can be overridden using the ``psf_stype_only`` option) and sources with SNRs that exceed ``psf_snrcut``. Next, a function is fit to determine how the size of sources (normalized by the median size) varies with the SNR. The function used is defined as :math:`\sigma / median = \sqrt(c_1^2 + c_2^2/SNR^2)`, where :math:`\sigma` is the size of the Gaussian and :math:`c_1` and :math:`c_2` are free parameters. Clipping of outliers is done during this fitting, controlled by the ``psf_nsig`` parameter. Lastly, unresolved sources are selected by choosing sources that lie within ``psf_kappa2`` times the rms of this best-fit sigma-SNR relation. As this last step can be unreliable for high-SNR sources, an additional selection can be made for the highest SNR sources using the ``psf_high_snr`` parameter. All sources with SNRs above ``psf_high_snr`` will be taken as unresolved.

* Next the image is tessellated using Voronoi tessellation to produce tiles within which the PSF shape is calculated (and assumed to be constant). The list of probable unresolved sources is filtered to select "calibrator" sources to use to determine the tessellation tiles. These sources are the brightest sources (known as the primary generators), defined as those sources that have SNRs in the top fraction of sources defined by ``psf_snrtop`` and that also have SNRs greater than ``psf_snrcutstack``. These sources are then grouped by their proximity, if they are within 50% of the distance to third closest source.

* The unresolved sources within each tile that have SNRs greater than ``psf_snrcutstack`` are then stacked to form a high-SNR PSF. For each tile, this PSF is fit with a Gaussian to recover its size. The significance of the variation in the sizes across the image is quantified.

* If the variation is significant, the major axis, minor axis, and position angle are then interpolated across the image. Where there is sufficient information, the interpolation is done using Delaunay triangulation; otherwise, the values within the tiles defined by tessellation are simply set to those of the appropriate PSF.

* Lastly, the deconvolved source sizes are adjusted to include the PSF variation as a function of position.

The options for this module are as follows:

.. parsed-literal::

    psf_vary_do ........... True : Calculate PSF variation across image
      :term:`psf_high_snr` ........ None : SNR above which all sources are taken to be
                                   unresolved. E.g., psf_high_snr = 20.0. None => no
                                   such selection is made
      :term:`psf_itess_method` ....... 0 : 0 = normal, 1 = 0 + round, 2 = LogSNR, 3 =
                                   SqrtLogSNR
      :term:`psf_kappa2` ........... 2.0 : Kappa for clipping for analytic fit
      :term:`psf_nsig` ............. 3.0 : Kappa for clipping within each bin
      :term:`psf_over` ............... 2 : Factor of nyquist sample for binning bmaj, etc. vs
                                   SNR
      :term:`psf_snrcut` .......... 10.0 : Minimum SNR for statistics
      :term:`psf_snrcutstack` ..... 15.0 : Unresolved sources with higher SNR taken for
                                   stacked psfs
      :term:`psf_snrtop` .......... 0.15 : Fraction of SNR > snrcut as primary generators
      :term:`psf_stype_only` ...... True : Restrict sources used in PSF variation
                                   estimating to be only of type 'S'

.. glossary::

    psf_high_snr
        This parameter is a float (default is ``None``). Gaussians with SNR greater than this are used to determine the PSF
        variation, even if they are deemed to be resolved. This corrects for the
        unreliability at high SNRs in the algorithm used to find unresolved
        sources. The minimum value is 20.0. If ``None``, then no such selection is made.

    psf_itess_method
        This parameter is an integer (default is 0) which can be 0, 1, 2 or 3, which
        corresponds to a tessellation method. If 0, 2 or 3, then the weights
        used for Voronoi tessellation are unity, log(SNR) and sqrt[log(SNR)]
        where SNR is the signal to noise ratio of the generator in a tile. If 1,
        then the image is tessellated such that each tile has smooth boundaries
        instead of straight lines, using pixel-dependent weights.

    psf_kappa2
        This parameter is a float (default is 2.0). When iteratively arriving at a statistically probable set of
        'unresolved' sources, the fitted major and minor axis sizes versus SNR
        are binned and fitted with analytical functions. Those Gaussians which
        are within ``psf_kappa2`` times the fitted rms from the fitted median are
        then considered 'unresolved' and are used further to estimate the PSFs.

    psf_nsig
        This parameter is a float (default is 3.0). When constructing a set of 'unresolved' sources for psf estimation, the
        (clipped) median, rms and mean of major and minor axis sizes of
        Gaussians versus SNR within each bin is calculated using ``kappa = psf_nsig``.

    psf_over
        This parameter is an integer (default is 2). When constructing a set of 'unresolved' sources for psf estimation, this parameter controls the factor of nyquist sample for binning bmaj, etc. vs SNR.

    psf_snrcut
        This parameter is a float (default is 10.0). Only Gaussians with SNR greater than this are considered for processing.
        The minimum value is 5.0

    psf_snrcutstack
        This parameter is a float (default is 15.0). Only Gaussians with SNR greater than this are used for estimating PSF
        images in each tile.

    psf_snrtop
        This parameter is a float (default is 0.15). If ``psf_generators`` is 'calibrator', then the peak pixels of Gaussians
        which are the ``psf_snrtop`` fraction of the SNR distribution are taken as Voronoi
        generators.

    psf_stype_only
        This parameter is a Boolean (default is ``False``). If ``True``\, sources are restricted to be only of type 'S'.

.. _spectralindex_do:

Spectral index module
---------------------
If ``spectralindex_do = True`` (and the input image has more than one frequency), then spectral indices are calculated for the sources in the following way:

* The rms maps for the remaining channels are determined.

* Neighboring channels are averages to attempt to obtain the target SNR per channel for a given source, set by the ``specind_snr`` parameter.

    .. note::

        No color corrections are applied during averaging. However, unless the source spectral index is very steep or the channels are very wide, the correction is minimal. See :ref:`colorcorrections` for details.

* Flux densities are measured for both individual Gaussians and for total sources. Once source flux densities have been measured in each channel, the SEDs are fit with a polynomial function. The best-fit parameters are then included in any catalogs that are written out (see :ref:`write_catalog`). In addition, plots of the fits can be viewed with the ``show_fit`` task (see :ref:`showfit`).

The options for this module are as follows:

.. parsed-literal::

    spectralindex_do ...... True : Calculate spectral indices (for multi-channel
                                   image)
      :term:`flagchan_rms` ........ True : Flag channels before (averaging and) extracting
                                   spectral index, if their rms if more than 5
                                   (clipped) sigma outside the median rms over all
                                   channels, but only if <= 10% of channels
      :term:`flagchan_snr` ........ True : Flag channels that do not meet SNR criterion set
                                   by specind_snr
      :term:`specind_maxchan` ........ 0 : Maximum number of channels to average for a
                                   given source when when attempting to meet target
                                   SNR. 1 => no averaging; 0 => no maximum
      :term:`specind_snr` .......... 3.0 : Target SNR to use when fitting power law. If
                                   there is insufficient SNR, neighboring channels
                                   are averaged to obtain the target SNR

.. glossary::

    flagchan_rms
        This parameter is a Boolean (default is ``True``). If ``True``, then the clipped rms and median (r and m) of the clipped rms of
        each channel is calculated. Those channels whose clipped rms is greater
        than 4r away from m are flagged prior to averaging and calculating
        spectral indices from the image cube. However, these channels are
        flagged only if the total number of these bad channels does not exceed
        10% of the total number of channels themselves.

    flagchan_snr
        This parameter is a Boolean (default is ``True``). If ``True``, then flux densities in channels that do not meet the target SNR are not used in fitting.

    specind_maxchan
        This parameter is an integer (default is 0) that sets the maximum number of channels that can be averaged together to attempt to reach the target SNR set by the ``specind_snr`` parameter. If 0, there is no limit to the number of channels that can be averaged. If 1, no averaging will be done.

    specind_snr
        This parameter is a float (default is 3.0) that sets the target SNR to use when fitting for the spectral index. If there is insufficient SNR, neighboring channels are averaged to obtain the target SNR. The maximum allowable number of channels to average is determined by the ``specind_maxchan`` parameter. Channels (after averaging) that fail to meet the target SNR are not used in fitting.

.. _polarisation_do:

Polarization module
-------------------
If ``polarisation_do = True``, then the polarization properties of the sources are calculated. First, if ``pi_fit = True``, source detection is performed on the polarized intensity (PI) image [#f3]_ to detect sources without Stokes I counterparts. The polarization module then calculates the I, Q, U, and V flux densities, the total, linear, and circular polarisation fractions and the linear polarisation angle of each Gaussian and source. The linear polarisation angle is defined from North, with positive angles towards East. Flux densities are calculated by fitting the normalization of the Gaussians found using the Stokes I or PI images.

For linearly polarised emission, the signal and noise add vectorially, giving a
Rice distribution instead of a Gaussian one. To correct for this, a bias
is estimated and removed from the polarisation fraction using the same method used for the
NVSS catalog (see ftp://ftp.cv.nrao.edu/pub/nvss/catalog.ps). Errors on the linear and total
polarisation fractions and polarisation angle are estimated using the debiased polarised flux density
and standard error propagation. See Sparks & Axon (1999) [#f4]_ for a more detailed treatment.

The options for this module are as follows:

.. parsed-literal::

    polarisation_do ....... True : Find polarisation properties
      :term:`pi_fit` .............. True : Check the polarized intesity (PI) image for
                                   sources not found in Stokes I
      :term:`pi_thresh_isl` ....... None : Threshold for PI island boundary in number
                                   of sigma above the mean. None => use thresh_isl
      :term:`pi_thresh_pix` ....... None : Source detection threshold for PI image:
                                   threshold for the island peak in number of sigma
                                   above the mean. None => use thresh_pix

.. glossary::

    pi_fit
        This parameter is a Boolean (default is ``True``). If ``True``, the polarized intensity image is searched for sources not
        present in the Stokes I image. If any such sources are found, they are
        added to the the Stokes I source lists. Use the ``pi_thresh_pix`` and
        ``pi_thresh_isl`` parameters to control island detection in the PI image.

    pi_thresh_isl
        This parameter is a float (default is ``None``) that determines the region to which fitting is done in the
        polarized intensity (PI) image. If ``None``, the value is set to that of the ``thresh_isl`` parameter. A higher value will produce smaller
        islands, and hence smaller regions that are considered in the fits. A
        lower value will produce larger islands. Use the ``pi_thresh_pix`` parameter
        to set the detection threshold for sources. Generally, ``pi_thresh_isl``
        should be lower than ``pi_thresh_pix``.

    pi_thresh_pix
        This parameter is a float (default is ``None``) that sets the overall detection threshold for islands in the
        polarized intensity (PI) image (i.e. pi_thresh_pix = 5 will find all
        sources with peak flux densities per beam of 5-sigma or greater). If ``None``, the value is set to that of the ``thresh_pix`` parameter. Use the ``pi_thresh_isl``
        parameter to control how much of each island is used in fitting.
        Generally, ``pi_thresh_pix`` should be larger than ``pi_thresh_isl``.

.. _shapelet_do:

Shapelet decomposition module
-----------------------------
If ``shapelet_do = True``, then islands are decomposed into shapelets. Shapelets are a set of 2-D basis functions (for details, see Refregier 2003 [#f5]_) that can be used to completely model any source, typically with far fewer parameters than pixels in the source. Shapelets are useful in particular for modeling complex islands that are not well modeled by Gaussians alone. PyBDSM can currently fit cartesian shapelets to an image. The shapelet parameters can be written to a catalog using ``write_catalog`` (see :ref:`write_catalog`).

For each island of emission, a shapelet decomposition is done after estimating the best values of the
center, the scale :math:`\beta`, and nmax in the following way. First, an initial guess of :math:`\beta` is taken as :math:`2\sqrt{[m2(x)m2(y)]}`,
where :math:`m2` is the second moment over the island, based on shapeelt analysis
of simulated images of resolved sources. Similarly, a guess for nmax is taken as the minimum
of 14, and maximum of 10 and :math:`2n + 2` where :math:`n=\sqrt{(n^2 + m^2)}/n_p^n - 1`, where (n, m) is the size of
the island and :math:`n^m_p` is the synthesized beam minor axis FWHM in pixels. This guess for nmax is
based partly on simulations and partly on the requirememts of computing time, number of
constraints, etc, for shapelet decomposition.

These initial values are then used to calculate the optimal central position around which
to decompose the island. First, for every pixel in the island, the coefficients c12 and c21
are computed assuming that pixel as the centre of expansion. Next, the zero crossings for
every vertical line of the c12 image and horizontal line of the c21 image are computed. The
intersection point of these two zero-crossing vectors is then taken as the proper centre of the
expansion for the image. If this procedure does not work, then the first moment is taken as
the center.

This updated center position is used to compute the optimal :math:`\beta`, which is taken as the value of
:math:`\beta` that minimises the residual rms in the island area. Using this :math:`\beta`, the center is computed
once more and the final shapelet deocmposition is then made.

The options for this module are as follows:

.. parsed-literal::

    shapelet_do ........... True : Decompose islands into shapelets
      :term:`shapelet_basis` .. 'cartesian': Basis set for shapelet decomposition:
                                   'cartesian' or 'polar'
      :term:`shapelet_fitmode` .... 'fit': Calculate shapelet coeff's by fitting ('fit') or
                                   integrating (None)

.. glossary::

    shapelet_basis
        This parameter is a string (default is ``'cartesian'``) that determines the type of shapelet
        basis used. Currently however, only cartesian is supported.

    shapelet_fitmode
        This parameter is a string (default is ``'fit'``) that determines the method of calculating
        shapelet coefficients. If ``None``, then these are calculated by integrating
        (actually, by summing over pixels, which introduces errors due to
        discretisation). If 'fit', then the coefficients are found by
        least-squares fitting of the shapelet basis functions to the image.

.. rubric:: Footnotes

.. [#f1] Condon, J. J. 1997, PASP, 109, 166

.. [#f2] Hopkins, A. M., Miller, C. J., Connolly, A. J., et al.  2002, AJ, 123, 1086

.. [#f3] The polarized intensity image is calculated as :math:`\sqrt{(Q^2 + U^2)}`.

.. [#f4] Sparks, W. B., & Axon, D. J. 1999, PASP, 111, 1298

.. [#f5] Refregier, A. 2003, MNRAS, 338, 35.
