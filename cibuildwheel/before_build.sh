#!/bin/bash -eux
#
# This script should be called by `cibuildwheel` in the `before-build` stage.
# Environment variables `BOOST_BUILD_DIR` and `BOOST_INSTALL_DIR` must be set.

# Install oldest supported numpy
function install_numpy
{
  pip install oldest-supported-numpy
}

# Build the Boost Python libraries
function build_boost_python
{
  nproc=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
  inc_dir=$(python -c 'import sysconfig as sc; print(sc.get_path("include"))')
  cd "${BOOST_BUILD_DIR}/boost"
  ./bootstrap.sh --prefix="${BOOST_INSTALL_DIR}" \
    --with-libraries=python
  ./b2 -j "${nproc}" \
    cxxflags="-fPIC -I${inc_dir}" \
    link=static,shared \
    warnings=off \
    install
}

set -o pipefail
# install_numpy
# build_boost_python

find /usr -name "libgfortran*" -ls 2>/dev/null