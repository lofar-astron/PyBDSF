**********************
Capabilities of PyBDSM
**********************

PyBDSM can be run on FITS images (using PyFITS [#f1]_) or CASA images (using pyrap [#f2]_), including 3-D and 4-D cubes, and can handle blanked image pixels. If a spectral cube is given, then all source extraction as well as other computation (psf variation, wavelet decomposition, etc.) are done on a collapsed 2-D stokes I image. Once sources have been identified, their spectral and polarisation properties are then extracted from the full cubes. If you need a full 3-D Gaussian decomposition, then DUCHAMP [#f3]_ is what you need.

PyBDSM performs the following tasks:

    * Reads in the image, collapses specific frequency channels, with weights, and produces a 'continuum' image (the 'ch0' image) for all polarisations. The Stokes I ch0 image is used for all further computation.
    
    * Preprocessing is done, whereby some basic parameters like image statistics are computed. Also, any input parameters that are left to default are calculated using sensible algorithms.
    
    * The background rms and mean images are computed. If the variation in these images is not statistically significant, then a constant value is taken. The parameters for this calculation are computed generically and hence do not have information about, for example, the typical size of the artifacts around bright sources. These parameters (e.g., :term:`rms_box`) are probably the only ones the user needs to take care to specify.
    
    * A constant threshold for separating source and noise pixels is set. This threshold can be either a hard threshold or calculated using the False Detection Rate algorithm.
    
    * Using these parameters, islands of contiguous source emission are identified. Islands are the basic units which are operated upon subsequently.
    
    * Each island is now fit with multiple Gaussians. Depending on the number of degrees of freedom, etc, the size of the fitted Gaussian could be fixed to be the restoring beam. The fitted Gaussians are then flagged to produce a list of acceptable set of Gaussians.
    
    * Each island can also be decomposed into shapelets. Currently only cartesian shapelets are implemented, and only one shapelet set, with the same scale in both dimensions, can be fit to an island. The shapelets parameters can be written out as ASCII or FITS tables.
    
    * Residual FITS images are computed, for both Gaussians and shapelets. The Gaussian parameters can be written out in various formats (ASCII, FITS tables, LOFAR BBS, ds9 region files, AIPS star, Kvis, etc). Shapelet parameters can be written out to FITS tables.
    
    * Gaussians within a given island are grouped into discrete sources.
    
    * If a frequency cube is input, then for each source identified in an island, the spectral index is computed. If possible, a spectral index is calculated for each Gaussian as well. This is done for point as well as extended sources.
    
    * If all four Stokes images are present, then the polarisation percentage and angle are calculated for each source. 
    
    * The residual ch0 image, after subtracting fitted Gaussians, is processed using the *à trous* wavelet transform to generate images at various scales. Islands are identified in each of these wavelet images and fitted with Gaussians, all of which are then grouped to form pyramidal sources. These can be used further by the user as a starting point for morphological filters.
    
    * Since the ionosphere affects low frequencies significantly, PyBDSM can also estimate the spatial variation of the PSF across the image, which can be used to correct various source parameters.

.. rubric:: Footnotes
.. [#f1] http://www.stsci.edu/resources/software_hardware/pyfits/
.. [#f2] http://code.google.com/p/pyrap/
.. [#f3] http://www.atnf.csiro.au/people/Matthew.Whiting/Duchamp/
