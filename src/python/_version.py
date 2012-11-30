"""Version module.

This module simply stores the version and svn revision numbers, as well
as a changelog. The svn revision number will be updated automatically
whenever there is a change to this file. However, if no change is made
to this file, the revision number will get out of sync. Therefore, one
must update this file with each (significant) update of the code:
adding to the changelog will naturally do this.
"""

# Version number
__version__ = '1.5.0'

# Store svn Revision number. For this to work, one also needs to do:
#
# "svn propset svn:keywords Revision CEP/PyBDSM/src/python/_version.py"
#
# from the LOFAR directory. Then, the revision number is
# added automatically with each update to this file. The line below does not
# need to be edited by hand.
__revision__ = filter(str.isdigit, "$Revision$")


# Change log
def changelog():
    """
    PyBDSM Changelog.
    -----------------------------------------------------------------------

    2012/11/30 - Fix to bypass bug in matplotlib when display variable
        is not set.

    2012/11/21 - Fixed bug that caused a crash when a detection image
        was used. Fixed a bug with incorrect save directory when
        "plot_allgaus" is True.

    2012/10/29 - Version 1.5.0

    2012/10/29 - Improved handling of WCS information so that a much
        greater variety of WCS systems may be used. Fixed a bug in logging
        that occurred when negative values were found in the rms map.
        Updated installation instructions.

    2012/10/12 - Version 1.4.5

    2012/10/12 - Added option ("incl_empty") to include empty islands (that
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
        0.0 Jy/beam unless "mean_map = 'const'" was explicitly specified.
        Fixed a bug in Gaussian fitting that could cause an island to be
        skipped.

    2012/10/02 - Fixed a bug in the PSF vary module that resulted in
        incorrect PSF generators being used. Added an option to smooth
        the resulting PSF images ("psf_smooth"). Parallelized the PSF
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

    2012/09/19 - Added "broadcast" option to show_fit task to send
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
        The number of cores to be used can be specified with the "ncores"
        option (default is to use all). Fixed bug in SED plotting in
        the show_fit task.

    2012/08/29 - Fixed incorrect terminal size in parameter listing. Added
        logging of non-default input parameters and internally derived
        parameters.

    2012/08/22 - Version 1.3.2

    2012/08/22 - Fixed a bug that caused the user-specified rms_box to be
        ignored. Added an option to enable the Monte Carlo error estimation
        for 'M'-type sources (the "do_mc_errors" option), which is now
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
        Aperture fluxes, if measured, are output in the "srl" catalogs.
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
    	minor axis maps in terms of sigma instead of FWHM. Added psf_vary
    	option (psf_stype_only) to allow PSF fitting to non- S-type sources
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
