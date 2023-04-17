#!/bin/bash -eux
#
# This script should be called by `cibuildwheel` in the `before-all` stage.
#
# This script will download and untar the Boost C++ source files.
# `${BOOST_VERSION}` should be set as `x.y.z`; the untarred sources can be
# found in the directory `${BOOST_BUILD_DIR}/boost`.
# Both environment variables must have been set.

# Download and untar the Boost C++ source files
# Rename the source directory to `boost`
function download_and_untar_boost
{
  major=$(echo "${BOOST_VERSION}" | cut -d. -f1)
  minor=$(echo "${BOOST_VERSION}" | cut -d. -f2)
  patch=$(echo "${BOOST_VERSION}" | cut -d. -f3)

  name="boost"
  long_name="${name}_${major}_${minor}_${patch}"
  site="https://sourceforge.net"
  directory="projects/${name}/files/${name}/${major}.${minor}.${patch}"
  file="${long_name}.tar.bz2"

  url="${site}/${directory}/${file}"

  rm -rf "${BOOST_BUILD_DIR}"
  mkdir -p "${BOOST_BUILD_DIR}"
  cd "${BOOST_BUILD_DIR}"
  curl -L -o - "${url}" | tar -xjf -
  mv "${long_name}" "${name}"
}

# find /usr -name "*gfortran*" -ls 2>/dev/null || true

set -o pipefail
download_and_untar_boost
