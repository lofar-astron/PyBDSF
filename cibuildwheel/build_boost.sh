#!/bin/bash -eux

PYMAJOR=$(python -c 'import sys; print(sys.version_info.major)')
PYMINOR=$(python -c 'import sys; print(sys.version_info.minor)')
PYUNICODE=$([ ${PYMAJOR} -eq 3 -a ${PYMINOR} -le 7 ] && echo "m" || echo "")
TARGET=cp${PYMAJOR}${PYMINOR}-cp${PYMAJOR}${PYMINOR}${PYUNICODE}
THREADS=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')

PY_INC_DIR=$(python -c 'import sysconfig as sc; print(sc.get_path("include"))')
PLATFORM=$(python -c 'import sys; print(sys.platform)')

echo "PLATFORM=${PLATFORM}"

# LDFLAGS="${LDFLAGS-} $(python-config --ldflags)"

which python
python --version

echo "BOOST_VERSION=${BOOST_VERSION-}"
echo "BOOST_BUILD_DIR=${BOOST_BUILD_DIR-}"
echo "BOOST_INSTALL_DIR=${BOOST_INSTALL_DIR-}"
echo "CFLAGS=${CFLAGS-}"
echo "LDFLAGS=${LDFLAGS-}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH-}"
echo "DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH-}"

# # find / -name "*fortran*" -ls 2>/dev/null
# exit 1

cd ${BOOST_BUILD_DIR}/boost
./bootstrap.sh --prefix=${BOOST_INSTALL_DIR} \
  --with-libraries=python
./b2 -j${THREADS} \
  cxxflags="-fPIC -I${PY_INC_DIR}" \
  link=static,shared install
# ./bootstrap.sh --prefix=${BOOST_INSTALL_DIR} \
#   --with-libraries=python \
#   --with-python=/opt/python/${TARGET}/bin/python \
#   --with-python-version=${PYMAJOR}.${PYMINOR} \
#   --with-python-root=/opt/python/${TARGET}
# ./b2 -j${THREADS} \
#   cxxflags="-fPIC -I/opt/python/${TARGET}/include/python${PYMAJOR}.${PYMINOR}${PYUNICODE}/" \
#   link=static,shared install
