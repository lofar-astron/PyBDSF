"""Changelog module.

This module records all the relevant changes made to the software.
The change list must be kept up-to-date manually.
"""

def changelog():
    """
    PyBDSF Changelog.
    -----------------------------------------------------------------------
    2023/05/22 - Version 1.10.3

    2023/05/08 - Fix build issue with Python 3.11 (#205)

    2023/05/03 - Use cibuildwheel to build binary wheels (#203)
        Build binary wheels for Linux and MacOS (Intel).
        Drop support for Python 3.6.

    2023/05/02 - Fix #198 (#199)
        Use the new method call `canvas.manager.set_window_title`

    2023/04/28 - Replace Travis CI with GitHub actions (#196)

    2023/02/10 - Version 1.10.2

    2023/02/10 - Fix issues with numpy versions >= 1.24 (#193)

    2022/11/28 - Switch to `manylinux2014` for building binary wheels (#191)

    2022/11/23 - Fix ImportError in setuptools (#190)

    2022/10/31 - Add binary wheels for Python 3.10 (#186)

    2022/10/14 - Fix various documentation issues (#185)

    2022/10/11 - Add logfilename option (#181)

    2022/10/05 - Use len() instead of numpy.alen() (#180)

    2022/02/14 - Version 1.10.1: Fix Numpy API incompatibility issue

    2022/02/09 - Version 1.10.0

    2022/02/09 - Update some functions as required by scipy versions >= 1.8.0
        (PR #172)

    2022/02/09 - Fix build issues with Python 3.8, end support for Python < 3.6,
        add support for Python 3.8 and 3.9, and make installation of the interactive
        pybdsf shell optional (PR #169)

    2022/02/09 - Improve handling of the beam in the spectral index module
        (PR #165)

    2021/05/05 - Improve handling of large, complex islands (PR #160)

    2020/04/07 - Allow a file to be supplied for the ch0 image used in the
        spectral index module (PR #127)

    2019/12/05 - Version 1.9.2

    2019/12/04 - Fix exception behaviour if spline order change does not work

    2019/09/27 - Add check for frequency info in header

    2019/09/25 - Version 1.9.1

    2019/09/25 - Fix various minor bugs

    2019/06/06 - Fix blank_limit check_low error (#100)

    2019/05/09 - Fix various shapelet decomposition issues

    2019/05/08 - Fix crash in Gaussian fitting (#96)

    2019/03/25 - Version 1.9.0

    2018/10/18 - Add support for Python 3

    2018/10/18 - Fix various minor bugs

    2018/10/12 - Version 1.8.15

    2018/10/09 - Fix segfault in Gaussian fitting (#63)

    2018/10/04 - Fix math domain error (#76)

    2018/06/21 - Fix setup.py for boost versions > 1.63

    2018/05/18 - Version 1.8.14

    2018/05/18 - Fix an error on total flux density (#50)

    2018/05/18 - Add the possibility to provide an external noise and mean maps (#43)

    2018/05/18 - Append the image FITS header into catalog FITS header (#53)

    2018/05/18 - Make PyBDSF compatible with newer boost libraries, specifically
        those used in Ubuntu 18.04 (#55)

    2017/11/17 - Version 1.8.13

    2017/11/17 - Remove deprecated boolean operators

    2017/09/01 - Version 1.8.12

    2017/09/01 - Fix crash with tiny regions

    2017/09/01 - Fix very low centroid peak fluxes

    2017/09/01 - Fix compile error with numpy 1.13

    2017/06/01 - Version 1.8.11

    2017/06/01 - Fix for interactive shell problem

    2017/05/31 - Version 1.8.10

    2017/05/31 - Fixes for various installation and runtime issues on modern systems.

    2017/03/23 - Version 1.8.9

    2017/03/23 - Fix to bug that causes an error when grouping Gaussians
        into sources

    2017/03/17 - Version 1.8.8

    2017/03/17 - Rename to PyBDSF, move to github, add setup.py installer

    2017/02/28 - Fix to issues related to numpy >= 1.12 and astropy >= 1.3

    2016/06/10 - Version 1.8.7

    2016/06/10 - Fix to bug that caused incorrect output images when input
        image was not square.

    2016/01/20 - Version 1.8.6

    2016/01/15 - Fix to bug that caused incorrect island mask when two
        islands are very close together.

    2015/12/07 - Fix to bug that caused crash when image is masked and
        the src_ra_dec option is used.

    2015/11/30 - Version 1.8.5

    2015/11/25 - Fix to bug in export_image that resulted in incorrect
        output image when both trim_box and pad_image were used.

    2015/11/20 - Fix to bug in wavelet module related to merging of islands.

    2015/11/20 - Fix to bug in polarization module related to numbering of
        new islands.

    2015/11/20 - Fix to bug in spectral index module related to rms map
        calculation.

    2015/11/20 - Added option to use much faster (but also much more memory
        intensive) SciPy fftconvolve function instead of custom PyBDSM one.
        The option (use_scipy_fft) defaults to True.

    2015/11/20 - Increased number of digits for values in output text
        catalogs

    2015/08/06 - Version 1.8.4

    2015/08/06 - Improved speed of wavelet module.

    2015/08/06 - Added option to use PyFFTW in wavelet module if available.

    2015/08/06 - Fix to IPython version check.

    2015/08/06 - Fix to bug that caused a failure to write shapelet models
        in FITS format.

    2014/11/07 - Fix to bug that caused a crash when both atrous_do = True
        and output_all = True. Fixed a bug that caused a crash on machines
        with only one core.

    2014/09/26 - Version 1.8.3

    2014/09/26 - Fix to bug that caused a crash when using the wavelet
        module and all Gaussians in an island were flagged.

    2014/07/03 - Mask will now be expanded to match input image shape. Fix
        to bug that caused image read failure when image lacks a Stokes axis.

    2014/05/14 - Version 1.8.2

    2014/05/15 - Fix to bug in CASA masks generated with export_image() that
        caused cleaning to fail in CASA 4.2 and above.

    2014/02/05 - Fix to bug that resulted in output file names being
        converted to lower case inappropriately.

    2014/01/14 - Version 1.8.1

    2014/01/13 - Added option (bbs_patches = 'mask') to allow patches in
        an output BBS sky model to be defined using a mask image.

    2014/01/09 - Fix to bug that caused the incl_empty option to be
        ignored when format='fits' in the write_catalog task.

    2013/12/05 - Enabled output of images in CASA format in the export_image
        task (img_format = 'casa'). Added an option to export_image task to
        export an island-mask image, with ones where there is emission and
        zeros elsewhere (image_type = 'island_mask'). Features in the island
        mask may be optionally dilated by specifying the number of dilation
        iterations with the mask_dilation parameter. Added an option to
        write a CASA region file to the write_catalog task (format =
        'casabox'). Added an option to write a CSV catalog to the
        write_catalog task (format = 'csv').

    2013/11/04 - Added error message when the rms is zero in some part of the
        rms map.

    2013/10/16 - Version 1.8.0

    2013/10/16 - Improved wavelet fitting. Added option so that wavelet
        fitting can be done to the sum of images on the remaining wavelet
        scales, improving the signal for fitting (controlled with the
        atrous_sum option). Added option so that user can choose whether to
        include new islands found only in the wavelet images in the final
        fit or not (controlled with the atrous_orig_isl option).

    2013/10/10 - Fixed a bug that could lead to incomplete fitting of
        some islands. Improved overall convergence of fits.

    2013/10/10 - Version 1.7.7

    2013/10/10 - Improved fitting of bright sources under certain
        circumstances.

    2013/09/27 - Version 1.7.6

    2013/09/27 - Changed caching behavior to ensure that temporary files
        are always deleted after they are no longer needed or on exit.

    2013/09/05 - Renamed blank_zeros to blank_limit. The blank_limit
        option now specifies a limit below which pixels are blanked.

    2013/09/05 - Enabled SAGECAL sky-model output.

    2013/09/02 - Version 1.7.5

    2013/09/02 - Fix to bug that caused a crash when images with 2 or
        3 axes were used. Improved rms and mean calculation (following the
        implementation used in PySE, see http://dare.uva.nl/document/174052
        for details). The threshold used to determine the clipped rms and
        mean values is now determined internally by default (kappa_clip =
        None).

    2013/08/27 - Version 1.7.4

    2013/08/29 - Fix to bug in show_fit() that caused error when
        'i' is pressed in the plot window and shapelet decomposition
        had not been done. Tweak to 'pybdsm' startup shell script to
        avoid problems with the Mac OS X matplotlib backend on non-
        framework Python installations (such as Anaconda Python).

    2013/08/28 - Fix to bug in process_image() that could result in
        wavelet Gaussians being excluded from model image under certain
        conditions.

    2013/08/27 - Version 1.7.3

    2013/08/27 - Fix to bug in image reading that caused images to be
        distorted.

    2013/08/23 - Version 1.7.2

    2013/08/23 - Improved handling of non-standard FITS CUNIT keywords.
        Improved loading of FITS images when trim_box is specified.

    2013/08/22 - Version 1.7.1

    2013/08/21 - Fix to bug that caused cached images to be deleted when
        rerunning an analysis. Fix to bug in show_fit() due to undefined
        images. Fix to bug in process_image() that would result in unneeded
        reprocessing.

    2013/08/20 - Version 1.7.0

    2013/08/19 - PyBDSM will now use Astropy if installed for FITS and WCS
        modules.

    2013/08/11 - Fix to avoid excessive memory usage in the wavelet module
        (replaced scipy.signal.fftconvolve with a custom function).

    2013/08/11 - Added option to use disk caching for internally derived
        images (do_cache). Caching can reduce memory usage and is
        therefore useful when processing large images.

    2013/07/11 - Implemented a variable operation chain for process_image
        (and img.process()) that allows unneeded steps to be skipped if
        the image is being reprocessed.

    2013/07/11 - Fixed a bug that could cause Gaussian fitting to hang
        during iterative fitting of large islands.

    2013/06/24 - Added option (fix_to_beam) to fix the size and position
        angle of Gaussians to the restoring beam during fitting. Fix to
        bug that caused the position angle used to initialize fitting to
        be incorrect.

    2013/03/22 - Version 1.6.1

    2013/03/21 - Fix to bug in ds9 and kvis catalog files that resulted in
        incorrect position angles. Fix to bug in position-dependent WCS
        transformations that caused incorrect source parameters in output
        catalogs. Added option to output uncorrected source parameters
        to a BBS sky model file (correct_proj).

    2013/03/14 - Removed sky transformations for flagged Gaussians, as
        these could sometimes give math domain errors. Disabled aperture
        flux measurement on wavelet images as it is not used/needed.

    2013/02/25 - Version 1.6.0

    2013/02/25 - Improved speed and accuracy of aperture flux
        calculation.

    2013/02/20 - Added option to use the curvature map method of
        Hancock et al. (2012) for the initial estimation of Gaussian
        parameters (ini_method = 'curvature') and for grouping of
        Gaussians into sources (group_method = 'curvature').

    2013/02/18 - Fix to bug in spectral index module that caused sources
        with multiple Gaussians to be skipped. Minor adjustments to the
        wavelet module to improve performance.

    2013/02/08 - Implemented position-dependent WCS transformations.

    2013/02/08 - Added option to fit to any arbitrary location in the
        image within a given radius (src_ra_dec and src_radius_pix).

    2013/02/04 - Fix to bug in wavelet module that caused crash when
       no Gaussians were fit to the main image.

    2013/01/30 - Fix to bug that resulted in incorrect numbering of
       wavelet Gaussians. Added 'srl' output in ds9 format when using
       output_all = True.

    2013/01/28 - Fix to bug in source grouping algorithm. Improved fitting
       when background mean is nonzero. Fix to allow images with GLAT and
       GLON WCS coordinates. Fix to bug when equinox is taken from the
       epoch keyword.

    2012/12/19 - Version 1.5.1

    2012/12/19 - Fix to bug in wavelet module that occurred when the
        center of the wavelet Gaussian lies outside of the image. Fix
        to re-enable srl output catalogs in ds9 region format. Fix to
        bug that resulted in the output directory not always being
        created. Added an option (aperture_posn), used when aperture
        fluxes are desired, to specify whether to center the aperture
        on the source centroid or the source peak.

    2012/12/02 - Changes to reduce memory usage, particularly in the
        wavelet module.

    2012/11/30 - Fix to bypass bug in matplotlib when display variable
        is not set.

    2012/11/21 - Fixed bug that caused a crash when a detection image
        was used. Fixed a bug with incorrect save directory when
        plot_allgaus = True.

    2012/10/29 - Version 1.5.0

    2012/10/29 - Improved handling of WCS information so that a much
        greater variety of WCS systems may be used. Fixed a bug in logging
        that occurred when negative values were found in the rms map.
        Updated installation instructions.

    2012/10/12 - Version 1.4.5

    2012/10/12 - Added option (incl_empty) to include empty islands (that
        have no un-flagged Gaussians) in output catalogs. Any such empty
        islands are given negative source IDs and positions given by the
        location of the peak of the island.

    2012/10/10 - Fixed a bug in Gaussian fitting that could cause a crash
        when fitting fails. Fixed a bug in parallelization that could
        cause a crash due to improper concatenation of result lists.

    2012/10/09 - Version 1.4.4

    2012/10/09 - Improved logging. Added a warning when one or more islands
        are not fit properly (i.e., no valid, unflagged Gaussians were
        fit). Fixed a bug in parallelization of Gaussian fitting that
        could cause a crash due to improper mapping of island lists to
        processes.

    2012/10/05 - Added code to handle images with no unblanked pixels.
        Improved fitting robustness.

    2012/10/04 - Version 1.4.3

    2012/10/04 - Fixed a bug in the mean map calculation that caused mean
        maps with constant values (i.e., non-2D maps) to have values of
        0.0 Jy/beam unless mean_map = 'const' was explicitly specified.
        Fixed a bug in Gaussian fitting that could cause an island to be
        skipped.

    2012/10/02 - Fixed a bug in the PSF vary module that resulted in
        incorrect PSF generators being used. Added an option to smooth
        the resulting PSF images (psf_smooth). Parallelized the PSF
        interpolation and smoothing steps. Improved PSF vary documentation.

    2012/09/25 - Version 1.4.2

    2012/09/25 - Dramatically reduced the time required to identify valid
        wavelet islands.

    2012/09/21 - Fixed bug that resulted in output FITS gaul tables being
        improperly sorted. Fixed cosmetic bug in the statusbar that could
        sometimes cause improper formatting. Added example of SAMP usage
        to the documentation.

    2012/09/20 - Version 1.4.1

    2012/09/20 - Fixed a bug in the wavelet module that caused a crash when
        no Gaussians were fit to the ch0 image.

    2012/09/19 - Added option (broadcast) to show_fit task to send
        coordinates and row highlight request to a SAMP hub when a Gaussian
        is clicked. Fixed bug in aperture flux masking that sometimes caused
        the mask to be the wrong shape.

    2012/09/18 - Added option to send images and catalogs to a SAMP hub
        (activated by setting outfile = 'SAMP' in the export_image and
        write_catalog tasks).

    2012/09/13 - Improved speed of plotting when images are large and in
        mean/rms map generation. Fixed bug that caused residual image
        statistics to fail when NaNs are present.

    2012/09/11 - Version 1.4.0

    2012/09/11 - Parallelized Gaussian fitting, shapelet decomposition,
        validation of wavelet islands, and mean/rms map generation.
        The number of cores to be used can be specified with the ncores
        option (default is to use up to 8). Fixed bug in SED plotting in
        the show_fit task.

    2012/08/29 - Fixed incorrect terminal size in parameter listing. Added
        logging of non-default input parameters and internally derived
        parameters.

    2012/08/22 - Version 1.3.2

    2012/08/22 - Fixed a bug that caused the user-specified rms_box to be
        ignored. Added an option to enable the Monte Carlo error estimation
        for 'M'-type sources (the do_mc_errors option), which is now
        disabled by default.

    2012/07/11 - Version 1.3.1

    2012/07/11 - Cleaned up unused options.

    2012/07/10 - Fixed a bug that caused a segfault during Gaussian
        fitting. Fixed a bug that caused a crash when a detection image
        is used.

    2012/07/05 - Fixed a bug that caused images written when output_all =
        True to be transposed. Added frequency information to all output
        images. Improved fitting robustness to prevent rare cases in
        which no valid Gaussians could be fit to an island. Modified the
        island-finding routine to handle NaNs properly.

    2012/07/03 - Version 1.3

    2012/07/03 - Fixed a bug in calculation of the positional errors of
        Gaussians. If interactive=True and image is large (> 4096 pixels),
        display is limited to 'ch0_islands' only; otherwise, show_fit()
        is very slow. Tweaked show_fit() to better display a single image.

    2012/07/02 - Adjusted rms_box algorithm to check for negative rms
        values (due to interpolation with cubic spline). If negative
        values are found, either the box size is increased or the
        interpolation is done with order=1 (bilinear) instead.

    2012/06/28 - Output now includes the residual image produced by
        using only wavelet Gaussians (if any) when atrous_do=True and
        output_all=True. Improved organization of files when
        output_all=True. Added logging of simple statistics (mean,
        std. dev, skew, and kurtosis) of the residual images.

    2012/06/22 - Included image rotation (if any) in beam definition.
        Rotation angle can vary across the image (defined by image WCS).

    2012/06/19 - Changed exception handling to raise exceptions when
        the interactive shell is not being used. Fixed bug that
        caused a crash when using show_fit() when no islands were
        found.

    2012/06/15 - Added Sagecal output format for Gaussian catalogs.

    2012/06/14 - Added check for newer versions of the PyBDSM
        software tar.gz file available on ftp.strw.leidenuniv.nl.

    2012/06/13 - Added total island flux (from sum of pixels) to
        "gaul" and "srl" catalogs.

    2012/06/06 - Version 1.2

    2012/06/06 - Added option to calculate fluxes within a specified
        aperture radius in pixels (set with the "aperture" option).
        Aperture fluxes, if measured, are output in the 'srl' catalogs.
        Changed code that determines terminal width to be more robust.

    2012/05/07 - Removed dependencies on matplotlib -- if matplotlib is
        not available, plotting is disabled. Corrected inconsistencies,
        spelling mistakes, etc. in help text and documentation. Cleaned
        up unneeded modules and files.

    2012/05/02 - Added option to output flux densities for every channel
        found by the spectral index module. Added option to spectral index
        module to allow use of flux densities that do not meet the desired
        SNR. Changed flag_maxsnr criterion to also flag if the peak flux
        density per beam of the Gaussian exceeds the value at its center.
        Removed incl_wavelet option.

    2012/04/20 - Promoted the adaptive_rms_box parameter to the main options
        listing and added the rms_box_bright option so that the user can
        specify either (or both) of the rms_boxes. Fixed bug in wavelet
        module so that invalid Gaussians (i.e., those that lie outside of
        islands in the ch0 image) are not used when making the residual
        images at each scale. Improved speed of Gaussian fitting to wavelet
        images. Fixed bug that resulted in pixels found to be outside the
        universe (check is enabled with the check_outsideuniv option) not
        being masked properly.

    2012/04/17 - Fixed bug in psf_vary module that resulted in PSF major and
        minor axis maps in terms of sigma instead of FWHM. Added option
        (psf_stype_only) to allow PSF fitting to non- S-type sources
        (useful if sources are very distorted).

    2012/04/12 - Fixed bug in adaptive scaling code that could cause
        incorrect small-scale rms_box size. Added a parameter
        (adaptive_thresh) that controls the minimum threshold for sources
        used to set the small-scale rms_box size.

    2012/04/02 - Implemented an adaptive scaling scheme for the rms_box
        parameter that shrinks the box size near bright sources and expands
        it far from them (enabled with the adaptive_rms_box option when
        rms_box=None). This scheme generally results in improved rms and
        mean maps when both strong artifacts and extended sources are
        present. Fixed bug that prevented plotting of results during wavelet
        decomposition when interactive = True.

    2012/03/29 - Fixed bug in wavelet module that could cause incorrect
        associations of Gaussians. Fixed bug in show_fit that displayed
        incorrect model and residual images when wavelets were used.

    2012/03/28 - Version 1.1

    2012/03/28 - Fixed bug that caused mask to be ignored when determining
        whether variations in rms and mean maps is significant. Fixed bug
        that caused internally derived rms_box value to be ignored.

    2012/03/27 - Modified calculation of rms_box parameter (when rms_box
        option is None) to work better with fields composed mainly of point
        sources when strong artifacts are present. Tweaked flagging on FWHM
        to prevent over-flagging of Gaussians in small islands. Changed
        wavelet module to flag Gaussians whose centers fall outside of
        islands found in the original image and removed atrous_orig_isl
        option (as redundant).

    2012/03/26 - Modified fitting of large islands to adopt an iterative
        fitting scheme that limits the number of Gaussians fit
        simultaneously per iteration to 10. This change speeds up fitting of
        large islands considerably. The options peak_fit and peak_maxsize
        control whether iterative fitting is done. Added new Gaussian
        flagging condition (flag_maxsize_fwhm) that flags Gaussians whose
        sigma contour times factor extends beyond the island boundary. This
        flag prevents fitting of Gaussians that extend far beyond the island
        boundary.

    2012/03/23 - Tweaked settings that affect fitting of Gaussians to
        improve fitting in general.

    2012/03/19 - Added output of shapelet parameters to FITS tables. Fixed
        issue with resizing of sources in spectral index module.

    2012/03/16 - Fixed bugs in polarisation module that caused incorrect
        polarization fractions.

    2012/03/13 - Improved plotting speed (by factor of ~ 4) in show_fit when
        there is a large number of islands. Simplified the spectral index
        module to make it more user friendly and stable. Added the option to
        use a "detection" image for island detection (the detection_image
        option); source properties are still measured from the main input
        image.

    2012/03/01 - Fixed a bug in the polarisation module that could result in
        incorrect flux densities. Changed logging module to suppress output
        of ANSI color codes to the log file.

    2012/02/27 - Implemented fitting of Gaussians in polarisation module,
        instead of simple summation of pixel values, to determine polarized
        flux densities.

    2012/02/17 - In scripts, process_image() will now accept a dictionary of
        parameters as input.

    2012/02/10 - Sources that appear only in Stokes Q or U (and hence not in
        Stokes I) are now identified and included in the polarisation
        module. This identification is done using the polarized intensity
        (PI) image. show_fit() and export_image() were updated to allow
        display and export of the PI image.

    2012/02/06 - Fixed bug in island splitting code that could result in
        duplicate Gaussians.

    2012/02/02 - Improved polarisation module. Polarization quantities are
        now calculated for Gaussians as well as sources.

    2012/01/27 - Fixed bug in psf_vary module that affected tesselation.
        Fixed many small typos in parameter descriptions.

    2012/01/18 - Fixed a bug that resulted in incorrect coordinates when the
        trim_box option was used with a CASA image. Added option
        (blank_zeros) to blank pixels in the input image that are exactly
        zero.

    2012/01/13 - Fixed minor bugs in the interactive shell and updated
        pybdsm.py to support IPython 0.12.

    2011/12/21 - Fixed bug in gaul2srl module due to rare cases in which an
        island has a negative rms value. Fixed a memory issue in which
        memory was not released after using show_fit.

    2011/11/28 - Added option to have minpix_isl estimated automatically as
        1/3 of the beam area. This estimate should help exclude false
        islands that are much smaller than the beam. This estimate is not
        let to fall below 6 pixels.

    2011/11/11 - Fixed bugs in source generation that would lead to masking
        of all pixels for certain sources during moment analysis. Adjusted
        calculation of jmax in wavelet module to use island sizes (instead
        of image size) if opts.atrous_orig_isl is True.

    2011/11/04 - Implemented new island fitting routine (enabled with the
        peak_fit option) that can speed up fitting of large islands. Changed
        plotting of Gaussians in show_fit to use Ellipse artists to improve
        plotting speed.

    2011/11/03 - Altered reading of images to correctly handle 4D cubes.
        Fixed bug in readimage that affected filenames.

    2011/10/26 - Extended psf_vary module to include fitting of stacked PSFs
        with Gaussians, interpolation of the resulting parameters across the
        image, and correction of the de- convolved source sizes using the
        interpolated PSFs. Changed plotting of Gaussians in show_fit() to
        use the FWHM instead of sigma. Modified error calculation of M
        sources to be more robust when sources are small. Fixed spelling of
        "gaussian" in bbs_patches option list.

    2011/10/24 - Many small bug fixes to the psf_vary module. Fixed use of
        input directory so that input files not in the current directory are
        handled correctly.

    2011/10/14 - Added residual rms and mean values to sources and source
        list catalogs. These values can be compared to background rms and
        mean values as a quick check of fit quality.

    2011/10/13 - Modified deconvolution to allow 1-D Gaussians and sources.
        Added FREQ0, EQUINOX, INIMAGE keywords to output fits catalogs.
        Fixed bug in source position angles. Adjusted column names of output
        catalogs slightly to be more descriptive.

    2011/10/12 - Added errors to source properties (using a Monte Carlo
        method for M sources). Fixed bug in output column names.

    2011/10/11 - Tweaked autocomplete to support IPython shell commands
        (e.g., "!more file.txt"). Fixed bug in gaul2srl that resulted in
        some very nearby Gaussians being placed into different sources.
        Added group_tol option so that user can adjust the tolerance of how
        Gaussians are grouped into sources.

    2011/10/05 - Added output of source lists. Changed name of write_gaul
        method to write_catalog (more general).

    2011/10/04 - Added option to force source grouping by island
        (group_by_isl). Added saving of parameters to a PyBDSM save file to
        Op_output.

    2011/09/21 - Fixed issue with shapelet centering failing: it now falls
        back to simple moment when this happens. Fixed issue with
        plotresults when shapelets are fit.

    2011/09/14 - Placed output column names and units in TC properties of
        Gaussians. This allows easy standardization of the column names and
        units.

    2011/09/13 - Fixes to trim_box and resetting of Image objects in
        interface.process(). Changed thr1 --> thr2 in fit_iter in
        guasfit.py, as bright sources are often "overfit" when using thr1,
        leading to large negative residuals. Restricted fitting of Gaussians
        to wavelet images to be only in islands found in the original image
        if opts.atrous_orig_isl is True.

    2011/09/08 - Version 1.0

    2011/09/08 - Versioning system changed to use _version.py.

    """
    pass
