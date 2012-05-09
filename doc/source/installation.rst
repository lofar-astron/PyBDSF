.. _installing:

**************************
Downloading and installing
**************************
.. note::

    If you are working on the LOFAR CEP I/II clusters, then PyBDSM is already installed. All that needs to be done is to initialize your environment as follows::
    
        $ use LofIm
        
Downloading the code
--------------------
The latest version of the code may be obtained from the LOFAR Subversion repository. First, change to the directory in which you wish to install PyBDSM. Then run the following command::

    $ svn co -N https://svn.astron.nl/LOFAR/trunk LOFAR
    $ svn up LOFAR/CMake    

Preparing to compile the code
-----------------------------
Before compiling the PyBDSM source code, you need to make sure you have the required dependencies:

* gfortran. Binaries are available from http://gcc.gnu.org/wiki/GFortranBinaries.
* Python (including numpy and scipy). The easiest way to install Python on a Mac is to use the 64-bit EPD Python distribution, available at http://enthought.com/products/epd.php.
* PyWCS. You can get PyWCS from https://trac6.assembla.com/astrolib.
* Boost. Get the latest version from www.boost.org. Only the Python libraries need to be compiled.

Compiling and installing
------------------------
Lastly, compile the software. To do so, change to the ``build/gnu_opt`` directory and execute ``make``::

    $ mkdir -p build/gnu_opt
    $ cd build/gnu_opt
    $ cmake -DBUILD_PACKAGES=PyBDSM -DUSE_LOG4CPLUS=OFF -DUSE_LOG4CXX=OFF ../../LOFAR
    $ make install
    
If successful, PyBDSM should now be installed in ``build/gnu_opt/installed/``. 

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
PyBDSM is currently under active development, and bug fixes and improvements are frequently committed to the Subversion repository. To update PyBDSM to the latest version, enter the following commands::

    $ cd $LOFAR/LOFAR
    $ svn update
    $ cd ../build/gnu_opt
    $ make install 
    
Major updates will be listed in :ref:`new`.
        

