"""Changelog module.

This module records all the relevant changes made to the software.
The change list must be kept up-to-date manually.
"""

def changelog():
    """
    PyBDSF Changelog.
    -----------------------------------------------------------------------

    2026/11/26 - Fix Monte Carlo uncertainty clipping for source positions and PAs (#585)

    2026/11/26 - Fix island peak clipping with varying noise (#590)

    2026/11/26 - Fix NaN propagation in `get_rot` (#592)

    2026/11/26 - Fix masked pixel handling when repairing negative RMS values (#577)

    2026/11/25 - Clean and optimize `psf_vary.py` (#552)

    2026/11/25 - Fix beam distortion when `fix_to_beam` is enabled (#570)

    2026/11/25 - Fix adaptive RMS message triggering (#593)

    2026/11/25 - Add info and warnings in the wavelet module (#597)

    2026/11/25 - Add a warning in `get_errors` when an exception occurs (#607)

    2026/11/25 - Add residual sigma statistics for standard and wavelet fitting (#600)

    2026/11/24 - Fix treshold warning triggering in `treshold.py` (#453)

    2026/11/24 - Add log message saying that source counts were extrapolated (#596)

    2026/11/24 - Fix automatic threshold selection for custom `thresh_pix` (#595)

    2026/11/23 - Add `residual_stats_do` to the docs (#602)

    2026/11/23 - Compute curvature map only once per island (#599)

    2026/11/23 - Fix DS9 generation modifies source parameters (#598)

    2026/11/23 - Add warning in `coords_to_isl` if user-supplied islands were skipped (#589)

    2026/11/23 - Fix silent external RMS map rejection (#580)

    2026/11/23 - Warn if fallback to 1.0 was executed in `fill_masked_regions` (#579)

    2026/11/23 - Add warning if adaptive RMS failed (#594)

    2026/11/23 - Update reference outputs after fixed Gaussian shape change (#604) (#606)

    2026/11/20 - Use Claude token

    2026/11/20 - Use fixed Gaussian shape for rank-deficient island masks (#604)

    2026/11/20 - Renamed reference model FITS files (#605)

    2026/11/20 - Add GitHub Actions workflow for @claude mentions (#582)

    2026/11/19 - Fix distorted beam detection for wrapped PAs (#571)

    2026/11/19 - Fix updating Gaussian counter after island merging (#490)

    2026/11/19 - Fix model filename (#484)

    2026/11/19 - Make `seterr` local only (#477)

    2026/11/17 - Fix peak flux comparison in curvature-based Gaussian grouping (#588)

    2026/11/17 - Fix Deprecation Warning in `gaul2srl` (#587)

    2026/11/17 - Fix normalization factor in the `norm_peak` calculation (#586)

    2026/11/16 - Fix incorrect peak-pixel coordinate calculation in `same_island_min` (#584)

    2026/11/16 - Fix `island_file` column header (#583)

    2026/11/16 - Fix order-dependent blending of small- and large-scale maps (#573)

    2026/11/11 - Fix coordinate mapping in RMS calculations (#561)

    2026/11/11 - Clean `gausfit.py` (#572)

    2026/11/11 - Merge pull request #569 from AlexKurek/Fix-FDR-crash-when-no-sources-are-found

    2026/11/11 - Merge pull request #544 from AlexKurek/Simplify-and-optimize-shapelets

    2026/11/09 - Combined reference updates (#576)

    2026/11/09 - Update reference outputs

    2026/11/09 - Update reference outputs and make output comparison stricter (#565)

    2026/11/09 - Fix wavelet output files (#548)

    2026/11/09 - Fix `inigaus_fbdsm` PA (#535)

    2026/11/09 - Fix intermittent hang of parallel_map when called from several threads (#566)

    2026/11/09 - Fix flaky test_interrupt_stops_workers on Python 3.10 (#567)

    2026/11/06 - Fix setting of `img.thresh_pix` when `pcrit == 0.0` (#559)

    2026/11/06 - Add microoptimization in `gausfit.py` (#558)

    2026/11/06 - Rewrite parallel_map with a process pool and dynamic scheduling (#503)

    2026/11/06 - Fix logging of compare outputs script (#563)

    2026/11/05 - Fix `fit_shapeletbasis` never running `leastsq` (#557)

    2026/11/05 - Fix message in polarization module (#556)

    2026/11/05 - Fix duplicated IDs for polarisation-only sources after island renumbering (#555)

    2026/11/05 - Skip Gaussian fit when adding a component fails (#554)

    2026/11/05 - Avoid unnecessary final Gaussian fit (#553)

    2026/11/05 - Add warning when automatic jmax is forced (#549)

    2026/11/04 - Took `m` out of the loop

    2026/11/04 - Fix filename completion for `tput` and `tget` (#547)

    2026/11/04 - Fix sorting in `output.py` (#546)

    2026/11/04 - Fix Gaussian numbering and preserve wavelet statistics (#542)

    2026/11/04 - Fix duplicated `gaus_num` after island renumbering in wavelet module (#541)

    2026/11/04 - Remove global warning ignoring from `plotresults.py` (#540)

    2026/11/04 - Remove 5 usnused methods and few variables (#539)

    2026/11/03 - Fix inconsistent Gaussian fitting thresholds (#538)

    2026/11/03 - Fix islands statistics after splitting (#537)

    2026/11/03 - Fix offset in `make_curvature_map` (#536)

    2026/11/02 - Modernize `readimage.py` (#534)

    2026/11/02 - Fix the 25-Gaussian limit in `inigaus_fbdsm` (#533)

    2026/11/02 - Remove C version of `bstat` (#525)

    2026/11/02 - Fix sorting of `p-values` in `threshold.py` (#530)

    2026/11/01 - Simplify and optimize shapelets

    2026/10/29 - Cache `rms_mean_map` results (#517)

    2026/10/29 - Fix masked pixels in nobeam Gaussian initialization (#500)

    2026/10/29 - Fix half-pixel offset in `generate_aperture` (#506)

    2026/10/29 - Clean output.py (#528)

    2026/10/29 - Fix masked neighbour handling in `get_maxima()` (#527)

    2026/10/29 - Remove Python 2 compatibility branches in `readimage.py` (#526)

    2026/10/27 - Fix incorrect axis coordinates in the padding edge case (#521)

    2026/10/27 - Fix detection of fully masked images in `rmsimage.py` (#520)

    2026/10/27 - Remove `morphfilter_pyramid` (#519)

    2026/10/26 - Speedup `rmsimage.py` (#514)

    2026/10/26 - Fix `outside_univ` (#513)

    2026/10/26 - Reduce RAM usage in wavelets (#512)

    2026/10/26 - Add early stop for `gaul` limit (#472)

    2026/10/26 - Fix `plist` passing to `get_errors` (#505)

    2026/10/26 - Fix inconsistent max position when the fit is rejected in `process_Multiple` (#502)

    2026/10/26 - Print PSF differences info (#501)

    2026/10/26 - Fix handling PA errors for circular Gaussians (#498)

    2026/10/26 - Fix extra islands criterion in `isl_tosplit` (#499)

    2026/10/26 - Fix background mean selection for polarization maps (#497)

    2026/10/26 - Fix background levels in Gaussian flagging (#496)

    2026/10/26 - Fix flux threshold consistency in `get_srcp` (#495)

    2026/10/26 - Revert "Fix automatic zero-mean mode handling in `rmsimage.py` (#493)" (#529)

    2026/10/23 - Fix automatic zero-mean mode handling in `rmsimage.py` (#493)

    2026/10/23 - Fix residual statistics for islands without Gaussians [bias] (#486)

    2026/10/22 - Fix updating of the Gaussian counter (#489)

    2026/10/22 - Add message in verbose fitting mode if the island was modified (#492)

    2026/10/22 - Fix island peak position for concave morphologies (#491)

    2026/10/22 - Fix the update of the island counter (#488)

    2026/10/22 - Fix source size errors assigment (#485)

    2026/10/21 - Remove output silencing from `functions.py` (#481)

    2026/10/21 - Remove `av` from Gaussian fitting (#480)

    2026/10/21 - Simplify gaussian fitting (#479)

    2026/10/21 - Call `momanalmask_gaus` only when needed (#478)

    2026/10/21 - Fix mean map test in `rmsimage.py` (#476)

    2026/10/21 - Remove dead code from `rmsimage.py` (#475)

    2026/10/21 - Speedup RMS calculations (#474)

    2026/10/21 - Remove unnecessary chceck in `deconv2` (#473)

    2026/10/21 - Fix eroding the island mask in `fit_island` (#465)

    2026/10/21 - Fix NaN support in `arrstatmask` (#424)

    2026/10/21 - Fix invalid IDs assigned to flagged Gaussians in `fgaul` (#412)

    2026/10/21 - Fix negative RMS map weights near image boundaries (#409)

    2026/10/20 - Fix `IndexError` and optimize bilinear interpolation in `fit_island` (#406)

    2026/10/20 - Speed up rms/mean map calculation in rmsimage (#504)

    2026/10/20 - Fix `ptbin` division in `psfvary.py` (#382)

    2026/10/20 - Add flagging of no-fit Gaussians (#483)

    2026/10/19 - Cache downloaded Boost tar-ball (#509)

    2026/10/18 - Cleanup scikit-build build directory (#508)

    2026/10/18 - Add NGC6048 FITS image via Git LFS (#507)

    2026/10/18 - Add benchmark script

    2026/10/12 - Fix WCS position error scaling for non-square pixels (#398)

    2026/10/09 - Skip `isl_tosplit` when `split_isl` is disabled (#471)

    2026/10/09 - Add warning on large islands before fitting (#467)

    2026/10/08 - Fix weighted sum of squares for the position angle in `psfvary.py` (#399)

    2026/10/08 - Fix underestimation of `deconv_size_skyE` (#397)

    2026/10/08 - Fix flux error aggregation for multi-component sources (#395)

    2026/10/08 - Fix determining the number of sources in `psfvary.py` (#374)

    2026/10/08 - Fix degree-to-pixel fallback in angdist2pixdist (#371)

    2026/10/05 - Speedup `inigaus_nobeam` (#464)

    2026/10/05 - Optimise Gaussian fitting (#463)

    2026/10/05 - Simplify `fit_mulgaus2d` (#462)

    2026/09/24 - Fix stale dof and mask in island fitting fallbacks (#461)

    2026/09/24 - Remove dead code in `_flag_gaussian` (#460)

    2026/09/24 - Remove `correct_borders` and `pad_array` (#457)

    2026/09/23 - Fix Gaussian deconvolution error handling in `deconv2` (#450)

    2026/09/22 - Replace `correct_borders` with `numpy.pad` (#455)

    2026/09/22 - Replace `pad_array` with `numpy.pad` (#454)

    2026/09/21 - Remove from __future__ import statements (#447)

    2026/09/17 - Fix `img.blankpix` calculation (#445)

    2026/09/17 - Fix flux and RMS calculations in plotting (#446)

    2026/09/15 - Fix support of peak position change in `rmsimage.py` (#419)

    2026/09/14 - Fix regression in `check_islands_for_overlap` (#440)

    2026/09/14 - Fix off by one error in island lookup for wavelet Gaussians in `wavelet_atrous.py` (#441)

    2026/09/11 - Revert `Fix RMS explosion near strong sources` (#438)

    2026/09/09 - Fix PA error calculation in `gaul2srl.py` (#429)

    2026/09/08 - Revert "Make p2s faster" (#434)

    2026/09/08 - Fix RMS explosion near strong sources (#431)

    2026/09/07 - Make progressbar more universal (#430)

    2026/09/07 - Fix bias towards 0 for flagged Gaussians (#425)

    2026/09/07 - Fix NaN support in `bstat` (#423)

    2026/09/07 - Add comment to `fit_island_iteratively` in `gausfit.py` (#427)

    2026/09/07 - Fix std calculation and zero division in moment (#422)

    2026/09/07 - Fix `s_imsize` dimension mismatch when sub-image reaches boundary (#421)

    2026/09/07 - Fix image boundary detection in `rms_mean_map` (#420)

    2026/09/07 - Fix masked region filling bias towards 0 and loop safety in `fill_masked_regions` (#415)

    2026/08/27 - Fix NaN propagation for very low SNR in `for_masked_mp` (#417)

    2026/08/27 - Fix negative `brightsize` values in `rmsimage.py` (#416)

    2026/08/27 - Fix NaN values handling in island peak detection (#414)

    2026/08/25 - Fix NaN contamination in island/source residual statistics `calc_resid_mean_rms` (#411)

    2026/08/24 - Fix shifting only along axis 0 in `imageshift` (#402)

    2026/08/24 - Fix `jmax` calcultion in wavelets (#381)

    2026/08/24 - Add another fallback in RMS MAD calculation (#376)

    2026/08/24 - Fix `area_pix` calculation in FDR (#358)

    2026/08/24 - Add warning to `treshold.py` (#346)

    2026/08/24 - Remove 14 unused functions from `functions.py` (#333)

    2026/08/24 - Fix threshold index search condition for ascending scflux (#408)

    2026/08/24 - Fix incorrect amplitude initialization in `fit_mask_1d` (#400)

    2026/08/21 - Fix stability of MC error estimation in `gaul2srl.py` (#396)

    2026/08/21 - Fix Gaussian size error scaling in `g2param_err` (#407)

    2026/08/21 - Remove unnecessary check in `isl_tosplit` (#405)

    2026/08/20 - Fix position angle normalization in `corrected_size` (#404)

    2026/08/20 - Add comment to `drawellipse` (#403)

    2026/08/20 - Fix crash due to typo in `fit_mulgaus2d` (#401)

    2026/08/20 - Fix Monte Carlo variance bias in `process_Multiple` (#394)

    2026/08/20 - Fix incorrect coordinate shift in `process_Multiple` (#393)

    2026/08/20 - Fix axis boundary checking in `same_island_min` (#392)

    2026/08/20 - Fix array index wrap-around in `same_island_min` (#391)

    2026/08/19 - Fix masked pixels included in image moment calculations (#389)

    2026/08/19 - Fix PA angle normalization in `deconv2` (#388)

    2026/08/19 - Replace `IOError` with `OSError` in `read_image_from_file` (#387)

    2026/08/19 - Fix mutation of `img.sources` causing duplicate entries (#386)

    2026/08/19 - Fix division in `gaul2srl.py` (#384)

    2026/08/19 - Fix `rms_box` size in wavelets (#383)

    2026/08/19 - Fix shapelets statistics (#380)

    2026/08/19 - Fix flux and RMS calculation in plotting (#379)

    2026/08/19 - Fix brightest pixel search for adaptive RMS map (#377)

    2026/07/27 - Fix typo in logging in `psfvary.py` (#375)

    2026/07/24 - Remove `nanmean` from functions.py (#297)

    2026/07/24 - Remove `fftconvolve` (#323)

    2026/07/24 - Fix island pixel counting in `rmsimage.py` validation check (#319)

    2026/07/24 - Clean `gaul2srl.py.py` (#318)

    2026/07/24 - Make `mapcoord_threaded` faster (#334)

    2026/07/24 - Make `p2s` faster (#359)

    2026/07/24 - Fix negative array dimensions in `psfvary.py` (#373)

    2026/07/24 - Fix plotting in `psfvary.py` (#372)

    2026/07/23 - Fix pixel count mismatch between total flux and its error (#370)

    2026/07/23 - Issue 366 (#368)

    2026/07/23 - Make `pixdist2angdist` faster (#282)

    2026/07/23 - Fix external mask in `list_and_sort_gaussians` (#331)

    2026/07/22 - Fix: mask pixels outside FOV before calculating statistics (#365)

    2026/07/22 - Clean `interface.py` (#348)

    2026/07/21 - Fix TOCTOU race in directory creation (#357)

    2026/07/21 - Fine-tune Gaussian fitting (#328)

    2026/07/21 - Remove duplicate from `opts.py` (#324)

    2026/07/21 - Fix `shapelet_getroot` logic (#321)

    2026/07/21 - Fix the calculation of mean island values (#351)

    2026/07/20 - Fix seconds formatting in `plotresults.py` (#364)

    2026/07/20 - Fix SED errorbars (#363)

    2026/07/20 - Fix `nmax` calculation in shapelets by removing incorrect upper limit (#362)

    2026/07/20 - Fix import in `shapefit.py` (#361)

    2026/07/20 - Implement mask support for `tresholdy.py` (#360)

    2026/07/17 - Fix formatting of `beam_spectrum` option in `print_opts()` (#349)

    2026/07/17 - Fix: subtract background when computing `total_flux` (#350)

    2026/07/17 - Add nproc() method (#354)

    2026/07/17 - Increase numeric tolerance (#355)

    2026/07/16 - Use correct mask

    2026/07/16 - Fix FDR algorithm (#345)

    2026/07/16 - Fix statistics of the residual image (#299)

    2026/07/16 - Remove `getTerminalSize` (#298)

    2026/07/16 - Fix island IDs map when `coords_to_isl` is used (#352)

    2026/07/16 - Fix `total_fluxE` computation (#301)

    2026/07/15 - Fix external mask in `list_and_sort_gaussians` (2) (#332)

    2026/07/15 - Improve the accuracy of `B1950toJ2000` (#330)

    2026/07/15 - Remove unnecessary `np.where` from `gausfit.py` (#329)

    2026/07/14 - Simplify adaptive RMS (#327)

    2026/07/14 - Fix `total_flux` computation (#295)

    2026/07/14 - Fix `RuntimeWarning: overflow encountered in cast` (#343)

    2026/07/13 - Fix initial Gaussian size estimation in `inigaus_nobeam` (#291)

    2026/07/13 - Fix computing the number of unmasked pixels assigned to current island in `Island` class (#294)

    2026/07/13 - Remove dead code from `treshold.py` (#344)

    2026/07/13 - Fix extended source filtering in adaptive rms box calculation (#339)

    2026/07/08 - Fix Numpy 2.5 deprecation (#341)

    2026/07/08 - Remove backport to `shutil_get_terminal_size` (#338)

    2026/07/08 - Simplify `for_masked` (#337)

    2026/07/08 - Fix Numpy 2.5 incompatibility (#325)

    2026/07/06 - Remove unused `np_where` from `rmsimage.py` (#315)

    2026/07/06 - More robust statistics of a very small RMS window (#335)

    2026/07/05 - Fine-tune inverse of the complementary error function in `treshold.py` (#336)

    2026/06/22 - Fix `try` block in `shape_findcen` (#322)

    2026/06/19 - Fix mean value computation in `rmsimage.py` (#320)

    2026/06/17 - Update reference output files (#314)

    2026/06/16 - Fix Gaussian fitting (#309)

    2026/06/16 - Clean `shapelets.py` (#312)

    2026/06/16 - Speedup statistics for the Gaussian residual image (#311)

    2026/06/15 - Remove duplicate from `gaul2srl.py` (#292)

    2026/06/15 - Cleanup of `rmsimage.py` (#289)

    2026/06/15 - Make printing of stats of residual images optional (#306)

    2026/06/15 - Don't use "input" symbol as it is a builtin Python function (#307)

    2026/06/15 - Remove `func_poly2d` and `func_poly2d_ini` (#285)

    2026/06/15 - Fix typo in `shapelet_check_centre` (#308)

    2026/06/12 - Fix indentation in `getzeroes_matrix` (#304)


    2026/06/12 - Version 1.14.1

    2026/06/12 - Use fork-style multiprocessing context (#305)


    2026/06/10 - Version 1.14.0

    2026/06/10 - Drop support for Python 3.9

    2026/06/09 - Remove dead code from `islands.py` (#293)

    2026/06/08 - Fix computing brightest pixel in the island (#296)

    2026/06/03 - Fix brightest source estimation (#290)

    2026/06/01 - Make statusbar faster (#286)

    2026/06/01 - Speedup `getTerminalSize` (#287)

    2026/05/29 - Remove duplicate in `opts` (#283)

    2026/05/29 - Make `check_islands_for_overlap` faster (#276)

    2026/05/22 - Fix renumbering of islands (#281)

    2026/05/21 - Remove `ntot_wvgaus` as a dead code (#280)

    2026/05/19 - Fix imports in `rmsimage.py` (#278)

    2026/05/18 - Remove unused imports from wavelet_atrous.py (#277)

    2026/05/01 - Make `gaus_2d` faster (#273)

    2026/04/09 - SCOOP-582 Optimize gausfit scheduling and Gaussian materialization (#271)

    2025/11/03 - Fix binary wheel build for Python 3.14 (#267)

    2025/10/30 - Fix test warnings (#266)

    2025/10/10 - Add RPATH

    2025/09/15 - Add argument types to function prototypes (#264)

    2025/09/02 - Merge pull request #260 from bnikolic/threadedmapcoordinates

    2025/09/01 - Merge pull request #262 from bnikolic/fastormasking

    2025/08/21 - Factor out selection of unmasked elements to improve perfomance

    2025/08/05 - Merge pull request #261 from bnikolic/stdevimprove

    2025/07/30 - Adapt also the last sigma calculation

    2025/07/30 - Save on standard deviation calculation

    2025/07/30 - pass the cores parameter

    2025/07/30 - Parallelise map_coordinates by multi-threaded map

    2025/07/23 - Improve performance of iterative sigma clipped rms (#259)

    2025/05/08 - Bump cibuildwheel version in GitHub workflow

    2025/05/06 - Check exit code of each joined process (#255)

    2025/04/10 - Fix failing publish to PyPI step in CI pipeline (#254)


    2025/04/08 - Version 1.13.0

    2025/04/08 - Issues/247 support numpy2.0 (#251)

    2025/04/07 - Fix exec() call for Python 3.13+ (#252)

    2024/11/25 - Version 1.12.0

    2024/11/19 - Fix string checks and clean up (#243)

    2024/11/15 - Fix unintended removal of handlers from the root logger (#245)

    2024/11/15 - Fix outdir being ignored (#244)

    2024/11/15 - rmsimage negative rms fallback, fixes #215 (#242)

    2024/10/11 - Fix bug that resulted in ncores being at most 8 and improve
                 method for getting no. of available cores (#240)

    2024/09/06 - Version 1.11.1

    2024/09/06 - Add macos arm binary wheels (#238)

    2024/09/04 - Remove manylinux2014 directory (#237)

    2024/09/04 - GitHubh action download artifact security update (#236)

    2024/07/26 - Version 1.11.0

    2024/07/26 - Prepare release v1.11 (#234)

    2024/07/01 - Add checks of output to tests (#231)

    2024/06/20 - Set upper limit on NumPy version (#230)

    2024/06/12 - Output island rms and mean values for wavelet Gaussians (#228)

    2024/05/06 - Use setuptools_scm (#224)

    2024/05/02 - Fix ReadtheDocs (#223)

    2024/04/24 - Remove pyfits version checks (#222)

    2024/03/21 - Fix issue with frequency_sp being ignored (#219)

    2024/03/11 - Fix python3.12 build (#217)

    2024/03/08 - Replace distutils for Python 3.12 (#216)

    2023/09/13 - Fix save of residual image and add option to save model image (#213)

    2023/07/13 - Fix issues with error estimation (#211)

    2023/07/05 - Fix various issue with the spectral index module (#210)

    2023/07/05 - Add options for use with a detection image (#209)

    2023/06/02 - Use scikit-build (#204)


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


    2022/02/14 - Version 1.10.1

    2022/02/14 - Fix Numpy API incompatibility issue


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
