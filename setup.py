#!/usr/bin/env python

from __future__ import print_function

import platform
from setuptools import Extension
from numpy.distutils.core import setup
import sys
import os
from ctypes.util import find_library
from os.path import join, realpath, dirname
import subprocess
import glob

def find_boost():
    # Find correct boost-python library.
    system = platform.system()
    if system == 'Linux':
        # Use version suffix if present
        boost_python = 'boost_python-py%s%s' % (sys.version_info[0], sys.version_info[1])
        if not find_library(boost_python):
            boost_python = "boost_python"
    elif system == 'Darwin':
        if sys.version_info[0] == 2:
            boost_python = "boost_python-mt"
        else:
            boost_python = "boost_python3-mt"
    return boost_python


boost_python = find_boost()

srcpath = join(dirname(realpath(__file__)),"src")

# Compile minpack, TODO: handle this with f2py
minpacksrc=["lmder.f", "lmpar.f", "qrfac.f", "qrsolv.f", "enorm.f", "dpmpar.f"]
for src in minpacksrc:
  subprocess.check_output(["gfortran", "-fPIC", "-c", src], cwd=join(srcpath, "minpack"))
subprocess.check_output("ar -cq libminpack.a *.o", cwd=join(srcpath, "minpack"), shell=True)

# Compile port3, TODO: handle this with f2py
port3src=["dnsg.f", "dn2g.f", "drnsg.f", "drn2g.f", "d1mach.f", "da7sst.f", "dc7vfn.f", "dd7tpr.f", "dd7upd.f", "df7hes.f", "dg7lit.f", "dg7qts.f", "ditsum.f", "divset.f", "dl7itv.f", "dl7ivm.f", "dl7mst.f", "dl7nvr.f", "dl7sqr.f", "dl7srt.f", "dl7svn.f", "dl7svx.f", "dl7tsq.f", "dl7tvm.f", "dl7vml.f", "dn2cvp.f", "dn2lrd.f", "dn2rdp.f", "do7prd.f", "dparck.f", "dq7apl.f", "dq7rad.f", "dq7rfh.f", "dr7mdc.f", "drldst.f", "ds7cpr.f", "ds7lup.f", "ds7lvm.f", "dv2axy.f", "dv2nrm.f", "dv7cpy.f", "dv7dfl.f", "dv7prm.f", "dv7scl.f", "dv7scp.f", "dv7swp.f", "i1mach.f", "i7mdcn.f", "stopx.f"]
for src in port3src:
  subprocess.check_output(["gfortran", "-fPIC", "-c", src], cwd=join(srcpath, "port3"))
subprocess.check_output("ar -cq libport3.a *.o", cwd=join(srcpath, "port3"), shell=True)

extensions=[]

fext=Extension(
    name="bdsf._pytesselate",
    sources=["src/fortran/pytess_simple.f",
             "src/fortran/pytess_roundness.f"]
    );
fext.f2py_options=[""]
extensions.append(fext)

extensions.append(Extension(
    name="bdsf._cbdsm",
    sources=["src/c++/Fitter_dn2g.cc",
             "src/c++/Fitter_dnsg.cc",
             "src/c++/Fitter_lmder.cc",
             "src/c++/MGFunction1.cc",
             "src/c++/MGFunction2.cc",
             "src/c++/cbdsm_main.cc",
             "src/c++/stat.cc",
             "src/c++/num_util/num_util.cpp"],
     libraries=['minpack', 'port3', 'gfortran', boost_python],
     include_dirs=["src/c++"],
     library_dirs=[join(srcpath,"minpack"), join(srcpath,"port3")]
     ))

extensions.append(Extension(
    name="bdsf.nat.natgridmodule",
    sources=glob.glob("natgrid/Src/*.c"),
    include_dirs = ["natgrid/Include"]
    ))

# HACK for supporting older versions of NumPy
for ext in extensions:
    ext.extra_f77_compile_args = []
    ext.extra_f90_compile_args = []

meta = dict(name='bdsf',
            version='1.8.12',
            author='David Rafferty',
            author_email='drafferty@hs.uni-hamburg.de',
            url='https://github.com/lofar-astron/PyBDSF',
            description='Blob Detection and Source Finder',
            long_description=open('README.rst', 'rt').read(),
            platforms='Linux, Mac OS X',
            packages=['bdsf', 'bdsf.nat'],
            package_dir = {'bdsf.nat': join('bdsf', 'nat')},
            classifiers=[
                'Intended Audience :: Science/Research',
                'Programming Language :: C++',
                'Programming Language :: Fortran',
                'Programming Language :: Python :: 2.7',
                'Topic :: Scientific/Engineering :: Astronomy'
                ],
            ext_modules=extensions,
            install_requires=['backports.shutil_get_terminal_size',
                              'numpy', 'scipy'],
            scripts = ['bdsf/pybdsf', 'bdsf/pybdsm'],
            zip_safe = False
            )

setup(**meta)
