.. _export_image:

**************************************************************
``export_image``: exporting internally derived images
**************************************************************

Internally derived images (e.g, the Gaussian model image) can be exported to FITS or CASA files using the ``export_image`` task:

.. parsed-literal::

    EXPORT_IMAGE: Write one or more images to a file.
    ================================================================================
    :term:`outfile` ............... None : Output file name. None => file is named
                                   automatically; 'SAMP' => send to SAMP hub (e.g., to
                                   TOPCAT, ds9, or Aladin)
    :term:`clobber` .............. False : Overwrite existing file?
    :term:`img_format` ........... 'fits': Format of output image: 'fits' or 'casa'
    :term:`img_type` ....... 'gaus_resid': Type of image to export: 'gaus_resid',
                                   'shap_resid', 'rms', 'mean', 'gaus_model',
                                   'shap_model', 'ch0', 'pi', 'psf_major', 'psf_minor',
                                   'psf_pa', 'psf_ratio', 'psf_ratio_aper', 'island_mask'
    :term:`mask_dilation` ............ 0 : Number of iterations to use for island-mask dilation.
                                   0 => no dilation
    :term:`pad_image` ............ False : Pad image (with zeros) to original size


Each of the parameters is described in detail below.

.. glossary::

    outfile
        This parameter is a string (default is ``None``) that sets the name of the output file. If ``None``, the file is named automatically. If 'SAMP' the image is sent to a running SAMP Hub (e.g., to ds9 or Aladin).

    clobber
        This parameter is a Boolean (default is ``False``) that determines whether existing files are overwritten or not.

    img_format
        This parameter is a string (default is ``'fits'``) that sets the output file format: ``'fits'`` - FITS format, ``'casa'`` - CASA format (requires pyrap).

    img_type
        This parameter is a string (default is ``'gaus_resid'``) that sets the type of image to export.
        The following images can be exported:

        * ``'ch0'`` - image used for source detection

        * ``'rms'`` - rms map image

        * ``'mean'`` - mean map image

        * ``'pi'`` - polarized intensity image

        * ``'gaus_resid'`` - Gaussian model residual image

        * ``'gaus_model'`` - Gaussian model image

        * ``'shap_resid'`` - Shapelet model residual image

        * ``'shap_model'`` - Shapelet model image

        * ``'psf_major'`` - image of major axis FWHM variation (arcsec)

        * ``'psf_minor'`` - image of minor axis FWHM variation (arcsec)

        * ``'psf_pa'`` - image of position angle variation (degrees east of north)

        * ``'psf_ratio'`` - image of peak-to-total flux variation (1/beam)

        * ``'psf_ratio_aper'`` - image of peak-to-aperture flux variation (1/beam)

        * ``'island_mask'`` - mask of islands (0 = outside island, 1 = inside island)

    mask_dilation
        This parameter is an integer (default is ``0``) that sets the number of dilation iterations to use when making the island mask. More iterations implies larger masked regions (one iteration expands the size of features in the mask by one pixel in all directions).

    pad_image
        This parameter is a Boolean (default is ``False``) that determines whether the output image is padded to be the same size as the original image (without any trimming defined by the ``trim_box`` parameter). If ``False``, the output image will have the size specified by the ``trim_box`` parameter.
