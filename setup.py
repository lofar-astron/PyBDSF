#!/usr/bin/env python

from __future__ import print_function

import platform
import setuptools
from numpy.distutils.core import setup, Extension
import numpy
import sys
from ctypes.util import find_library
from os.path import join, realpath, dirname
import subprocess
import glob
import distutils.cmd
import distutils.log
from numpy.distutils.command.build_ext import build_ext as numpy_build_ext
from distutils.command.clean import clean as distutils_clean
from distutils import ccompiler
import os
import warnings


no_boost_error = """
Could not find a Python boost library! Please use your package manager to install boost.

Or install it manually:

http://boostorg.github.io/python/doc/html/index.html
"""

minpack_src = ["lmder.f", "lmpar.f", "qrfac.f", "qrsolv.f", "enorm.f", "dpmpar.f"]
port3_src = ["dnsg.f", "dn2g.f", "drnsg.f", "drn2g.f", "d1mach.f", "da7sst.f",
             "dc7vfn.f", "dd7tpr.f", "dd7upd.f", "df7hes.f", "dg7lit.f", "dg7qts.f",
             "ditsum.f", "divset.f", "dl7itv.f", "dl7ivm.f", "dl7mst.f", "dl7nvr.f",
             "dl7sqr.f", "dl7srt.f", "dl7svn.f", "dl7svx.f", "dl7tsq.f", "dl7tvm.f",
             "dl7vml.f", "dn2cvp.f", "dn2lrd.f", "dn2rdp.f", "do7prd.f", "dparck.f",
             "dq7apl.f", "dq7rad.f", "dq7rfh.f", "dr7mdc.f", "drldst.f", "ds7cpr.f",
             "ds7lup.f", "ds7lvm.f", "dv2axy.f", "dv2nrm.f", "dv7cpy.f", "dv7dfl.f",
             "dv7prm.f", "dv7scl.f", "dv7scp.f", "dv7swp.f", "i1mach.f", "i7mdcn.f",
             "stopx.f"]
srcpath = join(dirname(realpath(__file__)), "src")


class CleanStatic(distutils.cmd.Command):
    """Custom command to remove static libraries and their intermediate files."""

    description = "remove static libs"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        for lib in ["minpack", "port3"]:
            self.announce(
                "Removing lib{}.a and object files".format(lib),
                level=distutils.log.INFO
            )
            subprocess.check_call(
                "rm -f *.o lib{}.a".format(lib),
                cwd=join(srcpath, lib),
                shell=True,
            )


class BuildStatic(CleanStatic):
    """Custom command to build static libraries."""

    # Compile minpack and port3, TODO: handle this with f2py
    description = "build static libs"

    def run(self):
        for lib, src_files in [["minpack", minpack_src], ["port3", port3_src]]:
            self.announce("Building lib{}".format(lib), level=distutils.log.INFO)
            for src in src_files:
                subprocess.check_output(
                    ["gfortran", "-fPIC", "-c", src],
                    cwd=join(srcpath, lib)
                )
            subprocess.check_output(
                "ar -cq lib{}.a *.o".format(lib),
                cwd=join(srcpath, lib),
                shell=True
            )


class BuildExt(numpy_build_ext):
    """Custom build_ext command that calls mbuild to build static libs"""

    def run(self):
        self.run_command("mbuild")
        numpy_build_ext.run(self)


class Clean(distutils_clean):
    """Custom clean command that calls mclean to clean static libs"""

    def run(self):
        self.run_command("mclean")
        distutils_clean.run(self)

def find_library_file(libname):
    """
    Try to get the directory of the specified library.
    It adds to the search path the library paths given to distutil's build_ext.
    """
    # Use a dummy argument parser to get user specified library dirs
    lib_dirs = [os.path.join(sys.prefix, 'lib'),
                             '/usr/local/lib',
                             '/usr/lib64',
                             '/usr/lib']
    try:
        lib_dirs.append(os.path.join('/usr/lib',
                                    getattr(sys, "implementation", sys)
                                    ._multiarch))
    except AttributeError:  # This is a non-multiarch aware Python.
        import sysconfig
        arch = sysconfig.get_config_var('MULTIARCH')
        if arch is not None:
            lib_dirs.append(os.path.join('/usr/lib', arch))

    if 'LD_LIBRARY_PATH' in os.environ:
        lib_dirs += os.environ['LD_LIBRARY_PATH'].split(':')

    compiler = ccompiler.new_compiler()
    return compiler.find_library_file(lib_dirs, libname)


def find_boost_python():
    """
    Find the name and path of boost-python

    Returns:
        library_name, e.g. 'boost_python-py36'        (a guess if boost is not found)
        library_dir,  e.g. '/opt/local/boost/lib'     ('' if boost is not found)
        include_dir,  e.g. '/opt/local/boost/include' ('' if boost is not found)
    """
    short_version = "{}{}".format(sys.version_info[0], sys.version_info[1])
    major_version = str(sys.version_info[0])

    # Prefer libraries with python version in their name over unversioned variants
    boostlibnames = ['boost_python-py' + short_version,
                     'boost_python' + short_version,
                     'boost_python' + major_version,
                     'boost_python',
                     ]
    # The -mt (multithread) extension is used on macOS but not Linux.
    # Look for it first to avoid ending up with a single-threaded version.
    boostlibnames = sum([[name + '-mt', name] for name in boostlibnames], [])

    for libboostname in boostlibnames:
        found_lib = find_library_file(libboostname)
        if found_lib:
            libdir = dirname(found_lib)
            includedir = join(dirname(libdir), "include")
            return libboostname, libdir, includedir

    warnings.warn(no_boost_error)
    return boostlibnames[0], '', ''


def find_boost_numpy():
    """
    Find the name and path of boost-numpy

    Returns:
        library_name, e.g. 'boost_numpy-py36'         (None if boost_numpy is not found)
        library_dir,  e.g. '/opt/local/boost/lib'     ('' if boost is not found)
        include_dir,  e.g. '/opt/local/boost/include' ('' if boost is not found)
    """
    short_version = "{}{}".format(sys.version_info[0], sys.version_info[1])
    major_version = str(sys.version_info[0])

    # Prefer libraries with python version in their name over unversioned variants
    boostlibnames = ['boost_numpy-py' + short_version,
                     'boost_numpy' + short_version,
                     'boost_numpy' + major_version,
                     'boost_numpy',
                     ]
    # The -mt (multithread) extension is used on macOS but not Linux.
    # Look for it first to avoid ending up with a single-threaded version.
    boostlibnames = sum([[name + '-mt', name] for name in boostlibnames], [])

    for libboostname in boostlibnames:
        found_lib = find_library_file(libboostname)
        if found_lib:
            libdir = dirname(found_lib)
            includedir = join(dirname(libdir), "include")
            return libboostname, libdir, includedir

    warnings.warn("No library boost_numpy found (this may be no problem)")
    return None, '', ''


def main():
    boost_python_libname, boost_python_libdir, boost_python_includedir = find_boost_python()
    boost_numpy_libname, boost_numpy_libdir, boost_numpy_includedir = find_boost_numpy()

    extensions = []

    fext = Extension(
        name="bdsf._pytesselate",
        sources=[
            "src/fortran/pytess_simple.f",
            "src/fortran/pytess_roundness.f"
        ]
    )
    fext.f2py_options = [""]
    extensions.append(fext)

    libraries = [
        'minpack', 'port3', 'gfortran', boost_python_libname
    ]
    if boost_numpy_libname is not None:
        libraries.append(boost_numpy_libname)

    extensions.append(Extension(
        name="bdsf._cbdsm",
        sources=[
            "src/c++/Fitter_dn2g.cc",
            "src/c++/Fitter_dnsg.cc",
            "src/c++/Fitter_lmder.cc",
            "src/c++/MGFunction1.cc",
            "src/c++/MGFunction2.cc",
            "src/c++/cbdsm_main.cc",
            "src/c++/stat.cc",
            "src/c++/num_util/num_util.cpp"
        ],
        libraries=libraries,
        include_dirs=[item for item in ("src/c++", boost_python_includedir, boost_numpy_includedir, numpy.get_include()) if item],
        library_dirs=[item for item in (join(srcpath, "minpack"), join(srcpath, "port3"), boost_python_libdir, boost_numpy_libdir) if item],
    ))

    extensions.append(Extension(
        name="bdsf.nat.natgridmodule",
        sources=glob.glob("natgrid/Src/*.c"),
        include_dirs=["natgrid/Include"]
    ))

    # HACK for supporting older versions of NumPy
    for ext in extensions:
        ext.extra_f77_compile_args = []
        ext.extra_f90_compile_args = []

    setup(
        name='bdsf',
        version='1.10.3a2',
        author='David Rafferty',
        author_email='drafferty@hs.uni-hamburg.de',
        url='https://github.com/lofar-astron/PyBDSF',
        description='Blob Detection and Source Finder',
        long_description=open('README.rst', 'rt').read(),
        long_description_content_type='text/x-rst',
        platforms='Linux, Mac OS X',
        packages=['bdsf', 'bdsf.nat'],
        package_dir={'bdsf.nat': join('bdsf', 'nat')},
        classifiers=[
            'Intended Audience :: Science/Research',
            'Programming Language :: C++',
            'Programming Language :: Fortran',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Topic :: Scientific/Engineering :: Astronomy'
        ],
        ext_modules=extensions,
        extras_require={
            'ishell': ['ipython<8.11', 'matplotlib']
        },
        install_requires=['backports.shutil_get_terminal_size',
                          'astropy', 'numpy', 'scipy'],
        entry_points = {
            'console_scripts': [
                'pybdsf = bdsf.pybdsf:main [ishell]',
                'pybdsm = bdsf.pybdsf:main [ishell]'
            ]
        },
        zip_safe=False,
        cmdclass={
            'mclean': CleanStatic,
            'mbuild': BuildStatic,
            'build_ext': BuildExt,
            'clean': Clean
        }
    )


if __name__ == "__main__":
    main()
