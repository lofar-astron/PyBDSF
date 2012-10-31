.. _showfit:

**************************************************
``show_fit``: visualizing the fit results
**************************************************

PyBDSM includes a task named ``show_fit`` that allows the user to quickly check the results of the ``process_image`` task. Use ``inp show_fit`` to list the parameters:

.. parsed-literal::

    SHOW_FIT: Show results of fit.
    ================================================================================
    :term:`broadcast` ............ False : Broadcast Gaussian and source IDs and coordinates
                                   to SAMP hub when a Gaussian is clicked?
    :term:`ch0_flagged` .......... False : Show the ch0 image with flagged Gaussians (if
                                   any) overplotted
    :term:`ch0_image` ............. True : Show the ch0 image. This is the image used for
                                   source detection
    :term:`ch0_islands` ........... True : Show the ch0 image with islands and Gaussians (if
                                   any) overplotted
    :term:`gmodel_image` .......... True : Show the Gaussian model image
    :term:`gresid_image` .......... True : Show the Gaussian residual image
    :term:`mean_image` ............ True : Show the background mean image
    :term:`pi_image` ............. False : Show the polarized intensity image
    :term:`psf_major` ............ False : Show the PSF major axis variation
    :term:`psf_minor` ............ False : Show the PSF minor axis variation
    :term:`psf_pa` ............... False : Show the PSF position angle variation
    :term:`rms_image` ............. True : Show the background rms image
    :term:`smodel_image` ......... False : Show the shapelet model image
    :term:`source_seds` .......... False : Plot the source SEDs and best-fit spectral
                                   indices (if image was processed with
                                   spectralindex_do = True). Sources may be chosen
                                   by ID with the 'c' key or, if ch0_islands = True,
                                   by picking a source with the mouse
    :term:`sresid_image` ......... False : Show the shapelet residual image

Each of the parameters is described in detail below.

.. glossary::

    broadcast
        This parameter is a Boolean (default is ``False``) that determines whether the Gaussian and source IDs and coordinates are sent to a running SAMP Hub when a Gaussian is clicked on. Note that for the IDs to be useful, a catalog must have been sent to the SAMP hub previously using the ``write_catalog`` task (with ``outfile = 'SAMP'``).

    ch0_flagged
        This parameter is a Boolean (default is ``False``) that determines whether to plot the ch0 image (the image used for source detection) with any flagged Gaussians overplotted.

    ch0_image
        This parameter is a Boolean (default is ``True``) that determines whether to plot the ch0 image (the image used for source detection).

    ch0_islands
        This parameter is a Boolean (default is ``True``) that determines whether to plot the ch0 image (the image used for source detection) with islands and Gaussians overplotted.

    gmodel_image
        This parameter is a Boolean (default is ``True``) that determines whether to plot the Gaussian model image.

    gresid_image
        This parameter is a Boolean (default is ``True``) that determines whether to plot the Gaussian residual image.

    mean_image
        This parameter is a Boolean (default is ``True``) that determines whether to plot the background mean image.

    pi_image
        This parameter is a Boolean (default is ``False``) that determines whether to plot the polarized intensity image.

    psf_major
        This parameter is a Boolean (default is ``False``) that determines whether to plot the variation of the major axis of the PSF.

    psf_minor
        This parameter is a Boolean (default is ``False``) that determines whether to plot the variation of the minor axis of the PSF.

    psf_pa
        This parameter is a Boolean (default is ``False``) that determines whether to plot the variation of the position angle of the PSF.

    rms_image
        This parameter is a Boolean (default is ``True``) that determines whether to plot the background rms image.

    smodel_image
        This parameter is a Boolean (default is ``False``) that determines whether to plot the shapelet model image.

    source_seds
        This parameter is a Boolean (default is ``False``) that determines whether to plot the source SEDs and best-fit spectral indices.

    sresid_image
        This parameter is a Boolean (default is ``False``) that determines whether to plot the shapelet residual image.
