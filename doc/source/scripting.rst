.. _scripting:

****************
Scripting PyBDSM
****************

Because PyBDSM is written in Python, it is straightforward to use PyBDSM non-interactively within Python scripts (for example, to automate source detection in a large number of images for which the optimal analysis parameters are known). To use PyBDSM in a Python script, import it by calling::

    import lofar.bdsm as bdsm

inside your script. 

.. note::

     If you are working on the LOFAR CEP I/II clusters, then you will also need to do::
    
        $ use Pythonlibs


Processing may then be done using ``process_image()`` as follows::

    img = bdsm.process_image(filename, <args>)          

where ``filename`` is the name of the image (in FITS or CASA format) or PyBDSM parameter save file and ``<args>`` is a comma-separated list of arguments defined as in the interactive environment (e.g., ``beam = (0.033, 0.033, 0.0), rms_map=False``). See :ref:`process_image` for details. 

.. note::

    The filename of the input image is also stored in the parameter save file. If you wish to override this filename (e.g., to use the saved parameters on a different image), give the save file as the first parameter and then explicitly set the filename. For example: ``img = bdsm.process_image('image1_savefile.sav', filename='image2.fits')``.

If the fit is successful, PyBDSM will return an Image object (in this example named ``img``) which contains the results of the fit (among many other things).  

When run in a Python script, it may be desirable to set ``output_all = True`` to write all output, including source lists, residual images, etc. to a directory named ``filename_pybdsm``. Optionally, the same tasks used in the interactive PyBDSM shell are available for examining the fit and writing out the source list, residual image, etc. These tasks are methods of the Image object returned by ``bdsm.process_image()`` and are described below. The input parameters to each of these tasks are the same as those available in the interactive shell (see the relevant task section for details).

``img.show_fit()``
    This method shows a quick summary of the fit by plotting the input image with the islands and Gaussians found, along with the model and residual images. See :ref:`showfit` for details.
    
``img.export_image()``
    Write an internally derived image (e.g., the model image) to a FITS file. See :ref:`export_image` for details.
    
``img.write_catalog()`` 
    This method writes the Gaussian or source list to a file. See :ref:`write_catalog` for details.

An example of using PyBDSM within a Python script is given in :ref:`script_example`. 
