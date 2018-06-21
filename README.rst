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
Installation is done through ``pip``, e.g. ``pip install .`` if you cloned the reposiory
or ``pip install https://github.com/lofar-astron/PyBDSF/archive/v1.8.15.tar.gz`` to install the latest release.

External requirements include the ubuntu packages (or similar packages in another Linux distribution):

* ``gfortran``
* ``libboost-python-dev``
* ``libboost-numpy-dev`` (Only if boost > 1.63)
* ``python-setuptools``.
  
Also, a working ``numpy`` installation is required.
At runtime, you will need ``scipy`` and either ``pyfits`` and ``pywcs`` or ``python-casacore`` or ``astropy``.

If you install as a user not using conda, use ``pip install --user``.
In this case, the script ``pybdsf`` is installed in ``~/.local/bin``, so you might want to add that to your ``$PATH``.

Installation on MacOS / OSX is more involved, you will need the packages mentioned above, for example installed with Homebrew.
You will need to tell `setup.py` to use the same compiler for fortran as for C++.

.. image:: https://travis-ci.org/lofar-astron/PyBDSF.svg?branch=master
    :target: https://travis-ci.org/lofar-astron/PyBDSF
