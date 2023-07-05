PyBDSF
======

PyBDSF (the Python **B**\ lob **D**\ etection and **S**\ ource **F**\ inder)
is a tool designed to decompose radio interferometry images into
sources and make available their properties for further use. PyBDSF can
decompose an image into a set of Gaussians, shapelets, or wavelets as
well as calculate spectral indices and polarization properties of
sources and measure the psf variation across an image. PyBDSF uses an
interactive environment based on CASA that will be familiar to most
radio astronomers. Additionally, PyBDSF may also be used in Python
scripts.

The documentation is currently hosted at https://pybdsf.readthedocs.io

Installation
------------
Installation can be done in a number of ways. In order of preference (read:
ease of use):

* Install the latest release from PyPI::

    pip install bdsf

  .. note:: The interactive shell ``pybdsf`` is no longer installed by default.
    To install it you have to specify the extra ``[ishell]``. For example::

      pip install bdsf[ishell]

* Install the ``master`` branch from the PyBDSF git repository::

    pip install git+https://github.com/lofar-astron/PyBDSF.git

  Or install a specific revision or release, for example ``v1.9.3``::

    pip install git+https://github.com/lofar-astron/PyBDSF.git@v1.9.3

* Install from a local source tree, e.g. after you cloned the git repository::

    pip install .

  or (to install the interactive shell as well)::

    pip install .[ishell]

If you get the error::

  RuntimeError: module compiled against API version 0xf but this version of numpy is 0xd

then please update ``numpy`` with ``pip install -U numpy``.

.. attention:: It is *not* recommend to use ``python setup.py install``. It is
  deprecated, and we do *not* support it.

External requirements include the ubuntu packages (or similar packages in another Linux distribution):

* ``gfortran``
* ``libboost-python-dev``
* ``libboost-numpy-dev`` (Only if boost > 1.63)
* ``python-setuptools``.

Also, a working ``numpy`` installation is required. At runtime, you will need ``scipy`` and either ``pyfits`` and ``pywcs`` or ``python-casacore`` or ``astropy``.

If you install as a user not using conda, use ``pip install --user``.
Make sure to use similar versions for gcc, g++ and gfortran
(use update-alternatives if multiple versions of gcc/g++/gfortran are present on the system).
In this case, the script ``pybdsf`` is installed in ``~/.local/bin``, so you might want to add that to your ``$PATH``.

Installation on MacOS / OSX is more involved, you will need the packages mentioned above, for example installed with Homebrew.
You will need to tell `setup.py` to use the same compiler for fortran as for C++. In case of problems, see https://github.com/lofar-astron/PyBDSF/issues/104#issuecomment-509267088 for some possible steps to try.

.. image:: https://github.com/lofar-astron/PyBDSF/actions/workflows/ci.yml/badge.svg?branch=master
    :target: https://github.com/lofar-astron/PyBDSF/actions/workflows/ci.yml
