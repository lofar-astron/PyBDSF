#!/bin/bash -eux
#
# This script should be called by `cibuildwheel` in the `before-build` stage.
#
# This script will first install the latest `numpy` 1.x version (2.x is not
# supported yet). Next the Boost Python libraries will be built from source,
# including the bindings to NumPy. The Boost sources must be in the directory
# `${BOOST_BUILD_DIR}/boost`. The libraries will be installed in the directory
# `${BOOST_INSTALL_DIR}`. Both environment variables must have been set.

# Ensure we start with a clean slate
function cleanup
{
  rm -rf "${BOOST_BUILD_DIR}/boost/bin.v2" 
  rm -rf "${BOOST_INSTALL_DIR}"
}

# Install latest numpy 1.x; we do not yet support numpy 2.x
function install_numpy
{
  pip install 'numpy<2'
}

# Build the Boost Python libraries
function build_boost_python
{
  nproc=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
  inc_dir=$(python -c 'import sysconfig as sc; print(sc.get_path("include"))')
  cd "${BOOST_BUILD_DIR}/boost"
  ./bootstrap.sh \
    --prefix="${BOOST_INSTALL_DIR}" \
    --with-toolset=gcc \
    --with-libraries=python
  ./b2 -j"${nproc}" -d0 \
    cxxflags="-fPIC -I${inc_dir}" \
    link=static,shared \
    install
}

set -o pipefail
cleanup
install_numpy
build_boost_python
