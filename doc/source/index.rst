====================
PyBDSF Documentation
====================

PyBDSF (the **Py**\thon **B**\lob **D**\etector and **S**\ource **F**\inder) is a tool designed to decompose radio interferometry images into sources and make available their properties for further use. PyBDSF can decompose an image into a set of Gaussians, shapelets, or wavelets as well as calculate spectral indices and polarization properties of sources and measure the psf variation across an image. PyBDSF uses an interactive environment based on CASA [#f1]_ that will be familiar to most radio astronomers. Additionally, PyBDSF may also be used in Python scripts.


.. toctree::
   :caption: Introduction
   :maxdepth: 2

   context
   capabilities


.. toctree::
   :caption: Obtaining PyBDSF
   :maxdepth: 2

   installation
   whats_new


.. toctree::
   :caption: User's Guide
   :maxdepth: 3

   ug_basics
   process_image
   show_fit
   export_image
   write_catalog
   scripting
   parameters


.. toctree::
   :caption: Analysis Examples
   :maxdepth: 2

   examples


.. toctree::
   :caption: Details of the Algorithms
   :maxdepth: 2

   algorithms


.. rubric:: Footnotes
.. [#f1] http://casa.nrao.edu
