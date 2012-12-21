.. _new:

**********
What's New
**********

Version 1.5.1 (2012/12/19):

    * Fix to bug in wavelet module that occurred when the center of the wavelet Gaussian lies outside of the image.

    * Fix to re-enable srl output catalogs in ds9 region format.

    * Fix to bug that resulted in the output directory not always being created.

    * Added an option (``aperture_posn``), used when aperture fluxes are desired, to specify whether to center the aperture on the source centroid or the source peak.

    * Changes to reduce memory usage, particularly in the wavelet module.

    * Fix to bypass bug in matplotlib when display variable is not set.

    * Fixed bug that caused a crash when a detection image was used.

    * Fixed a bug with incorrect save directory when "plot_allgaus" is True.

Version 1.5.0 (2012/10/29):

    * Improved WCS handling. PyBDSM can now read images with a much greater variety of WCS systems (e.g., the ``VOPT`` spectral system).

    * Fixed a bug related to the use of a detection image when a subimage is specified (with ``trim_box``).

Version 1.4.5 (2012/10/12):

    * Added option (``incl_empty``) to include empty islands (that have no un-flagged Gaussians) in output catalogs. Any such empty islands are given negative source IDs and have positions given by the location of the peak of the island.

    * Fixed a bug in Gaussian fitting that could cause a crash when fitting fails.

    * Fixed a bug in parallelization that could cause a crash due to improper concatenation of result lists.

Version 1.4.4 (2012/10/09):

    * Fixed a bug related to the parallelization of Gaussian fitting that could cause a crash due to improper mapping of island lists to processes.

    * Improved logging.

    * Added a warning when one or more islands are not fit (i.e., no valid, unflagged Gaussians were found).

    * Added code to handle images with no unblanked pixels.

    * Improved fitting robustness.

Version 1.4.3 (2012/10/04):

    * Fixed a bug in the mean map calculation that caused mean maps with constant values (i.e., non-2D maps) to have values of 0.0 Jy/beam unless ``mean_map = 'const'`` was explicitly specified.

    * Fixed a bug in the PSF vary module that resulted in incorrect PSF generators being used. Added an option to smooth the resulting PSF images (``psf_smooth``). Parallelized the PSF interpolation and smoothing steps. Improved PSF vary documentation.

Version 1.4.2 (2012/09/25):

    * Dramatically reduced time required to identify valid wavelet islands. Fixed bug that resulted in output FITS gaul tables being improperly sorted.

Version 1.4.1 (2012/09/11):

    * Added SAMP (Simple Application Messaging Protocol) support to the write_catalog, export_image, and show_fit tasks. These tasks can now use SAMP to communicate with other programs connected to a SAMP hub (e.g., ds9, Topcat, Aladin).

Version 1.4.0 (2012/09/11):

    * Parallelized Gaussian fitting, shapelet decomposition, validation of wavelet islands, and mean/rms map generation. The number of cores to be used can be specified with the ``ncores`` option (default is to use all).

Version 1.3.2 (2012/08/22):

    * Fixed a bug that could cause the user-specified ``rms_box`` value to be ignored. Added an option to enable the Monte Carlo error estimation for 'M'-type sources (the ``do_mc_errors`` option), which is now disabled by default.

Version 1.3.1 (2012/07/11):

    * Fixed a bug that caused images written when ``output_all = True`` to be transposed. Added frequency information to all output images. Improved fitting robustness to prevent rare cases in which no valid Gaussians could be fit to an island. Modified the island-finding routine to handle NaNs properly.

Version 1.3.0 (2012/07/03):

    * Fixed a bug in the calculation of positional errors for Gaussians.

    * Adjusted ``rms_box`` algorithm to check for negative rms values (due to interpolation with cubic spline). If negative values are found, either the box size is increased or the interpolation is done with ``order=1`` (bilinear) instead.

    * Output now includes the residual image produced using only wavelet Gaussians (if any) when ``atrous_do=True`` and ``output_all=True``.

    * Improved organization of files when ``output_all=True``.

    * Added logging of simple statistics (mean, std. dev, skew, and kurtosis) of the residual images.

    * Included image rotation (if any) in beam definition. Rotation angle can vary across the image (defined by image WCS).

    * Added Sagecal output format for Gaussian catalogs.

    * Added check for newer versions of the PyBDSM software ``tar.gz`` file available on ftp.strw.leidenuniv.nl.

    * Added total island flux (from sum of pixels) to ``gaul`` and ``srl`` catalogs.

Version 1.2 (2012/06/06):

    * Added option to output flux densities for every channel found by the spectral index module.

    * Added option to spectral index module to allow use of flux densities that do not meet the desired SNR.

    * Implemented an adaptive scaling scheme for the ``rms_box`` parameter that shrinks the box size near bright sources and expands it far from them (enabled with the ``adaptive_rms_box`` option when ``rms_box`` is None). This scheme generally results in improved rms and mean maps when both strong artifacts and extended sources are present.

    * Improved speed of Gaussian fitting to wavelet images.

    * Added option to calculate fluxes within a specified aperture radius in pixels (set with the ``aperture`` option). Aperture fluxes, if measured, are output in the ``srl`` format catalogs.

Version 1.1 (2012/03/28):

    * Tweaked settings that affect fitting of Gaussians to improve fitting in general.

    * Modified calculation of the ``rms_box`` parameter (when the ``rms_box`` option is None) to work better with fields composed mainly of point sources when strong artifacts are present.

    * Modified fitting of large islands to adopt an iterative fitting scheme that limits the number of Gaussians fit simultaneously per iteration to 10. This change speeds up fitting of large islands considerably.

    * Added the option to use a "detection" image for island detection (the ``detection_image`` option); source properties are still measured from the main input image. This option is particularly useful with images made with LOFAR's AWImager, as the uncorrected, flat-noise image (the ``*.restored`` image) is better for source detection than the corrected image (the ``*.restored.corr`` image).

    * Modified the polarization module so that sources that appear only in Stokes Q or U (and hence not in Stokes I) are now identified. This identification is done using the polarized intensity (PI) image.

    * Improved the plotting speed (by a factor of many) in ``show_fit`` when there are a large number of islands present.

    * Simplified the spectral index module to make it more user friendly and stable.

    * Altered reading of images to correctly handle 4D cubes.

    * Extended the ``psf_vary`` module to include fitting of stacked PSFs with Gaussians, interpolation of the resulting parameters across the image, and correction of the deconvolved source sizes using the interpolated PSFs.

    * Added residual rms and mean values to source catalogs. These values can be compared to background rms and mean values as a quick check of fit quality.

    * Added output of shapelet parameters as FITS tables.

    * Fixed many minor bugs.

See the changelog (accessible from the interactive shell using ``help changelog``) for details of all changes since the last version.
