PyBDSF
======

PyBDSF (the Python **B**\ lob **D**\ etection and **S**\ ource **F**\ inder, formerly
PyBDSM) is a tool designed to decompose radio interferometry images into
sources and make available their properties for further use. PyBDSF can
decompose an image into a set of Gaussians, shapelets, or wavelets as
well as calculate spectral indices and polarization properties of
sources and measure the psf variation across an image. PyBDSF uses an
interactive environment based on CASA that will be familiar to most
radio astronomers. Additionally, PyBDSF may also be used in Python
scripts.

The documentation is currently hosted at http://www.astron.nl/citt/pybdsf

Installation
------------
Installation is done through ``python setup.py``. External requirements include the ubuntu packages ``gfortran``, ``libboost-python-dev``, ``python-setuptools`` (or similar packages in another Linux distribution). Also, a working ``numpy`` installation is required.
