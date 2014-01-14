.. _installing:

**************************
Downloading and installing
**************************
.. note::

    If you are working on the LOFAR CEP I/II clusters, then PyBDSM is already installed. All that needs to be done is to initialize your environment as follows::

        $ use LofIm

Downloading the code
--------------------
The latest version of the code may be obtained as a gzipped tar file from the STRW FTP server at ftp://ftp.strw.leidenuniv.nl/pub/rafferty/PyBDSM (e.g., ``PyBDSM-1.8.0.tar.gz``). Once downloaded, extract the files in the directory where you would like to install PyBDSM. The files are all contained in a subdirectory named ``LOFAR``.

Preparing to compile the code
-----------------------------
Before compiling the PyBDSM source code, you need to make sure you have the required dependencies:

.. note::

    The minimal set of dependencies is usually part of most Linux distributions, so if you are installing PyBDSM on Linux you can likely skip to the next step (compiling and installing). On a Mac, you will also need to have XCode installed (from the Mac App Store), including the command-line tools (installed either from XCode's Preferences or, on 10.9 Mavericks, by running ``xcode-select --install`` in a terminal).

* Python 2.6 or 2.7 (including NumPy, SciPy, Matplotlib, and IPython). The easiest way to install Python and all of the required modules is to use the free 64-bit Anaconda distribution, available at http://www.continuum.io/downloads (Anaconda also includes Astropy, which is needed by PyBDSM). Python 3 is not yet supported.
* gfortran. Binaries are available from http://gcc.gnu.org/wiki/GFortranBinaries.
* A C++ compiler. Note that the default system compiler on OS 10.9 does not work with PyBDSM at this time, so it is necessary to install a recent version of the GCC compiler suite (e.g., the GCC 4.8 binaries from http://hpc.sourceforge.net). The easiest way to use these alternative compilers is to replace the system versions of the compilers in /usr/bin/ (i.e., cc, gcc, g++, c++) with these versions (before compiling Boost).
* Boost. Get the latest version from http://www.boost.org. Only the Python libraries need to be compiled. For example, on a Mac, do the following (which assumes the latest version is ``boost_1_49_0.tar.gz``)::

    $ cd /usr/local/
    $ sudo tar --bzip2 -xf ~/Downloads/boost_1_49_0.tar.gz
    $ cd boost_1_49_0/
    $ sudo ./bootstrap.sh --with-libraries=python
    $ sudo ./b2 install


.. note::

    If you don't have superuser access, you can install Boost to a local directory by adding::

        --prefix=path/to/installation/prefix

    to the bootstrap.sh command above and then passing this directory to the cmake command below by adding::

        -DBOOST_ROOT_DIR=/path/to/boost


* Astropy (if you use the Anaconda Python distribution above, Astropy is already included). You can get Astropy from http://www.astropy.org.


Compiling and installing
------------------------
Lastly, compile the software. To do so, change to the ``LOFAR`` directory and make a ``build/gnu_opt`` directory, go there, and execute ``make``::

    $ cd LOFAR
    $ mkdir -p build/gnu_opt
    $ cd build/gnu_opt
    $ cmake -DBUILD_PACKAGES=PyBDSM -DUSE_LOG4CPLUS=OFF -DUSE_LOG4CXX=OFF ../..
    $ make install

If successful, PyBDSM should now be installed in ``LOFAR/build/gnu_opt/installed/``.

.. _add_to_path:

Adding PyBDSM to your PATH
--------------------------
You can add PyBDSM to your PATH by adding the following lines to your ``.cshrc`` (for the C-shell) or ``.bash_profile`` files (for the Bash shell):

For the C-shell::

    setenv LOFAR <root directory of code tree>
    source $LOFAR/build/gnu_opt/installed/lofarinit.csh

For the Bash shell::

    export LOFAR="<root directory of code tree>"
    source $LOFAR/build/gnu_opt/installed/lofarinit.sh

.. note::

     If you are working on the LOFAR CEP I/II clusters, then you need only to do::

        $ use LofIm

Keeping up-to-date
------------------
PyBDSM is currently under active development, and bug fixes and improvements are frequently implemented. PyBDSM will automatically check for updates each time the interactive shell is started. To update PyBDSM to the latest version, download the new version and repeat the "compiling and installing" steps.

Major updates will be listed in :ref:`new`.


