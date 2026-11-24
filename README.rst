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

.. attention:: It is *not* recommend to use ``python setup.py install``. It is
  deprecated, and we do *not* support it.

PyBDSF is a pure Python package: no compilers are needed. Its
dependencies (``numpy``, ``scipy`` and ``astropy``) are installed by ``pip``.

If you install as a user not using conda, use ``pip install --user``.
In this case, the script ``pybdsf`` is installed in ``~/.local/bin``, so you might want to add that to your ``$PATH``.

.. image:: https://github.com/lofar-astron/PyBDSF/actions/workflows/ci.yml/badge.svg?branch=master
    :target: https://github.com/lofar-astron/PyBDSF/actions/workflows/ci.yml
