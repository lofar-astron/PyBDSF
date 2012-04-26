.. PyBDSM documentation master file, created by
   sphinx-quickstart on Thu Jan 19 13:27:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================
PyBDSM Documentation
====================

PyBDSM (the **Py**\thon **B**\lob **D**\etection and **S**\ource **M**\easurement software) is a tool designed to decompose radio interferometry images into sources and make available their properties for further use. PyBDSM can decompose an image into a set of Gaussians, shapelets, or wavelets as well as calculate spectral indices and polarization properties of sources and measure the psf variation across an image. PyBDSM uses an interactive environment based on CASA [#f1]_ that will be familiar to most radio astronomers. Additionally, PyBDSM may also be used in Python scripts.

.. .. image:: overview_image.png
..    :align: center
    

Introduction
============

.. toctree::
   :maxdepth: 2

   context
   capabilities


Obtaining PyBDSM
================

.. toctree::
   :maxdepth: 2

   installation
   whats_new


User's Guide
============

.. toctree::
   :maxdepth: 2

   ug_basics
   process_image
   show_fit
   export_image
   write_catalog
   scripting
   parameters
   

Analysis Examples
=================

.. toctree::
   :maxdepth: 2

   examples


Details of the Algorithms
=========================

.. toctree::
   :maxdepth: 2

   algorithms
   
.. rubric:: Footnotes
.. [#f1] http://casa.nrao.edu