.. _export_image:

**************************************************************
``export_image``: exporting internally derived images
**************************************************************

Internally derived images (e.g, the Gaussian model image) can be exported to FITS files using the ``export_image`` task:

.. parsed-literal::

    EXPORT_IMAGE: Write one or more images to a file.
    ================================================================================
    :term:`outfile` ............... None : Output file name. None => file is named
                                   automatically; 'SAMP' => send to SAMP hub (e.g., to
                                   TOPCAT, ds9, or Aladin)
    :term:`clobber` .............. False : Overwrite existing file?
    :term:`img_format` ........... 'fits': Format of output image: 'fits' or 'casa' (at the
                                   moment only 'fits' is supported)
    :term:`img_type` ....... 'gaus_resid': Type of image to export: 'gaus_resid',
                                   'shap_resid', 'rms', 'mean', 'gaus_model',
                                   'shap_model', 'ch0', 'pi'

Each of the parameters is described in detail below.

.. glossary::

    outfile
        This parameter is a string (default is ``None``) that sets the name of the output file. If ``None``, the file is named automatically. If 'SAMP' the image is sent to a running SAMP Hub (e.g., to ds9 or Aladin).

    clobber
        This parameter is a Boolean (default is ``False``) that determines whether existing files are overwritten or not.

    img_format
        This parameter is a string (default is ``'fits'``) that sets the output file format: ``'fits'`` - FITS format, ``'casa'`` - CASA format.

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

