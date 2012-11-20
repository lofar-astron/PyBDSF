.. _installing:

**************************
Downloading and installing
**************************
.. note::

    If you are working on the LOFAR CEP I/II clusters, then PyBDSM is already installed. All that needs to be done is to initialize your environment as follows::

        $ use LofIm

Downloading the code
--------------------
The latest version of the code may be obtained as a gzipped tar file from the STRW FTP server at ftp://ftp.strw.leidenuniv.nl/pub/rafferty/PyBDSM (e.g., ``PyBDSM-1.2.tar.gz``). Once downloaded, extract the files in the directory where you would like to install PyBDSM. The files are all contained in a subdirectory named ``LOFAR``.

Preparing to compile the code
-----------------------------
Before compiling the PyBDSM source code, you need to make sure you have the required dependencies:

* Python 2.6 or newer (including NumPy, SciPy, Matplotlib, and IPython). The easiest way to install Python and all of the required modules is to use the 64-bit EPD Python distribution, available at http://enthought.com/products/epd.php. For academic users, a free version is available at http://www.enthought.com/products/edudownload.php.
* gfortran. Binaries are available from http://gcc.gnu.org/wiki/GFortranBinaries.
* PyWCS. You can get PyWCS from https://trac6.assembla.com/astrolib.
* Boost. Get the latest version from http://www.boost.org. Only the Python libraries need to be compiled. For example, on a Mac, do the following (which assumes the latest version is ``boost_1_49_0.tar.gz``)::

    $ cd /usr/local/
    $ sudo tar --bzip2 -xf ~/Downloads/boost_1_49_0.tar.gz
    $ cd boost_1_49_0/
    $ sudo ./bootstrap.sh --with-libraries=python
    $ sudo ./b2 install


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
PyBDSM is currently under active development, and bug fixes and improvements are frequently implemented. PyBDSM will automatically check for updates each time the interactive shell is started. To update PyBDSM to the latest version, download the new version and repeat the above steps.

Major updates will be listed in :ref:`new`.


