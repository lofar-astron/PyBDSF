#!/bin/bash -eux
#
# This script should be called by `cibuildwheel` in the `before-build` stage.
#
# This script will first install the oldest supported `numpy` to maximize
# portability. Next the Boost Python libraries will be built from source,
# including the bindings to NumPy. The Boost sources must be in the directory
# `${BOOST_BUILD_DIR}/boost`. The libraries will be installed in the directory
# `${BOOST_INSTALL_DIR}`. Both environment variables must have been set.

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
    --with-libraries=python \
    --with-toolset=gcc
  ./b2 -j"${nproc}" \
    cxxflags="-fPIC -I${inc_dir}" \
    link=static,shared \
    install
}

# find /usr -name "libgfortran*" -ls 2>/dev/null || true

echo "CC=${CC}"
echo "CXX=${CXX}"
echo "FC=${FC}"

${CC} --version
${CXX} --version
${FC} --version

set -o pipefail
install_numpy
build_boost_python
