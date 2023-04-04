#!/bin/bash -eux

PYMAJOR=$(python -c 'import sys; print(sys.version_info.major)')
PYMINOR=$(python -c 'import sys; print(sys.version_info.minor)')
PYUNICODE=$([ ${PYMAJOR} -eq 3 -a ${PYMINOR} -le 7 ] && echo "m" || echo "")
TARGET=cp${PYMAJOR}${PYMINOR}-cp${PYMAJOR}${PYMINOR}${PYUNICODE}
THREADS=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')

# rm -rf /build
# mkdir /build
cd ${BOOST_BUILD_DIR}/boost
# tar xjf /src/boost.tar.bz2
./bootstrap.sh --prefix=${BOOST_INSTALL_DIR} \
  --with-libraries=python \
  --with-python=/opt/python/${TARGET}/bin/python \
  --with-python-version=${PYMAJOR}.${PYMINOR} \
  --with-python-root=/opt/python/${TARGET}
ls -l
./b2 -j${THREADS} \
  cxxflags="-fPIC -I/opt/python/${TARGET}/include/python${PYMAJOR}.${PYMINOR}${PYUNICODE}/" \
  link=static,shared install
