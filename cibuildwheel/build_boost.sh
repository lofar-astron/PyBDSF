#!/bin/bash -ex

PYMAJOR=$(python -c 'import sys; print(sys.version_info.major)')
PYMINOR=$(python -c 'import sys; print(sys.version_info.minor)')
PYUNICODE=$([ ${PYMAJOR} -eq 3 -a ${PYMINOR} -le 7 ] && echo "m" || echo "")
TARGET=cp${PYMAJOR}${PYMINOR}-cp${PYMAJOR}${PYMINOR}${PYUNICODE}
THREADS=$(nproc)

rm -rf /opt/boost /build
mkdir -p /src /build
cd /src
[ -f boost.tar.bz2 ] || \
    curl https://ufpr.dl.sourceforge.net/project/boost/boost/1.76.0/boost_1_76_0.tar.bz2 --output boost.tar.bz2
cd /build
tar jxf /src/boost.tar.bz2
cd boost_1_76_0
./bootstrap.sh --prefix=/opt/boost \
  --with-libraries=python \
  --with-python=/opt/python/${TARGET}/bin/python \
  --with-python-version=${PYMAJOR}.${PYMINOR} \
  --with-python-root=/opt/python/${TARGET}
./b2 -j${THREADS} \
  cxxflags="-fPIC -I/opt/python/${TARGET}/include/python${PYMAJOR}.${PYMINOR}${PYUNICODE}/" \
  link=static,shared install
