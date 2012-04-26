.. _new:

**********
What's New
**********

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