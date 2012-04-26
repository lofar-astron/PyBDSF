.. _installing:

**************************
Downloading and installing
**************************
.. note::

    If you are working on the LOFAR CEP I/II clusters, then PyBDSM is already installed. All that needs to be done is to initialize your environment as follows::
    
        $ use LUS
        $ use LofIm
        
    The order of the above commands is important, as path variable are not set correctly if the order is reversed.
    

Downloading the code
--------------------
The latest version of the code may be obtained from the LOFAR User's Software Group (USG) Subversion repository. First, change to the directory in which you wish to install PyBDSM. Then run the following command::

    $ svn co http://usg.lofar.org/svn/code/trunk lofarsoft

This command will call Subversion to check out the source code of PyBDSM, which will be placed in the lofarsoft directory. After the download is done, change to this directory::

    $ cd lofarsoft


Preparing to compile the code
-----------------------------
Before compiling the PyBDSM source code, you need to set up your environment correctly. To do so, enter the following commands:

For the C-shell::

    $ setenv LOFARSOFT <root directory of code tree>
    $ source $LOFARSOFT/devel_common/scripts/init.csh

For the Bash shell::

    $ export LOFARSOFT="<root directory of code tree>"
    $ source $LOFARSOFT/devel_common/scripts/init.sh

Although the installation script will attempt to automatically build and install any missing software that PyBDSM depends on, there are some unresolved issues that affect installation on Mac OS X. Therefore, if you are using a Mac, you need to install the following software first:

* gfortran. Binaries are available from http://gcc.gnu.org/wiki/GFortranBinaries.
* Python (including numpy and scipy). The easiest way to install Python on a Mac is to use the 64-bit EPD Python distribution, available at http://enthought.com/products/epd.php.
* PyWCS. You can get PyWCS from https://trac6.assembla.com/astrolib.

Next, run the bootstrap configuration script that will scan you system to determine if the needed dependencies are available::

    $ ./bootstrap

Compiling and installing
------------------------
Lastly, compile the software. To do so, change to the ``build`` directory and execute ``make``::

    $ cd build
    $ make anaamika

.. note::

    Anaamika is the name of the larger package that includes PyBDSM (see http://www.strw.leidenuniv.nl/~mohan/anaamika for details), but is no longer begin maintained. By default, the ``make`` command will build only the PyBDSM part of Anaamika.

Depending on your system, you may need to force the use of gfortran as follows::

    $ make anaamika FC=gfortran

If successful, PyBDSM should now be installed in ``$LOFARSOFT/release``. 

.. _add_to_path:

Adding PyBDSM to your PATH
--------------------------
You can add this directory to your PATH by adding the following lines to your ``.cshrc`` (for the C-shell) or ``.bash_profile`` files (for the Bash shell):

For the C-shell::

    setenv LOFARSOFT <root directory of code tree>
    source $LOFARSOFT/devel_common/scripts/init.csh

For the Bash shell::

    export LOFARSOFT="<root directory of code tree>"
    source $LOFARSOFT/devel_common/scripts/init.sh
    
.. note::

     If you are working on the LOFAR CEP I/II clusters, then you need only to do::
    
        $ use LUS
        $ use LofIm


Keeping up-to-date
------------------
PyBDSM is currently under active development, and bug fixes and improvements are frequently committed to the Subversion repository. To update PyBDSM to the latest version, enter the following commands::

    $ cd $LOFARSOFT
    $ svn update
    $ cd build
    $ make anaamika
    
Major updates will be listed in :ref:`new`.
        

