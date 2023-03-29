#!/bin/bash -ex
#
# This script must be run by `cibuildwheel`, before the actual build starts.
# See also `pyproject.toml`, section `[tool.cibuildwheel]`

yum install -y gcc-gfortran # boost-python36-devel

# download and build boost
# env 1>&2
# which python 1>&2
# python --version 1>&2
# ls -l /project 1>&2
# exit 1

# mkdir build
# cd build
# curl https://ufpr.dl.sourceforge.net/project/boost/boost/1.76.0/boost_1_76_0.tar.bz2 --output boost.tar.bz2
# tar jxf boost.tar.bz2
# cd boost_1_76_0
# ./bootstrap.sh --prefix=/opt/boost \
#   --with-libraries=python \
#   --with-python=/opt/python/${TARGET}/bin/python \
#   --with-python-version=${PYMAJOR}.${PYMINOR} \
#   --with-python-root=/opt/python/${TARGET}
# ./b2 -j${THREADS} \
#   cxxflags="-fPIC -I/opt/python/${TARGET}/include/python${PYMAJOR}.${PYMINOR}${PYUNICODE}/" \
#   link=static,shared install