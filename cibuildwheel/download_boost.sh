#!/bin/bash -eux

BOOST_VERSION="1.76.0"

BOOST_MAJOR=$(echo ${BOOST_VERSION} | cut -d. -f1)
BOOST_MINOR=$(echo ${BOOST_VERSION} | cut -d. -f2)
BOOST_PATCH=$(echo ${BOOST_VERSION} | cut -d. -f3)

NAME="boost"
LONG_NAME="${NAME}_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}"
BASE_URL="https://sourceforge.net"
DIRECTORY="projects/${NAME}/files/${NAME}/${BOOST_MAJOR}.${BOOST_MINOR}.${BOOST_PATCH}"
FILE="${LONG_NAME}.tar.bz2"

URL="${BASE_URL}/${DIRECTORY}/${FILE}"

rm -rf ${BOOST_BUILD_DIR}
mkdir -p ${BOOST_BUILD_DIR}
cd ${BOOST_BUILD_DIR}
curl -L -o - ${URL} | tar -xjf -
mv ${LONG_NAME} ${NAME}
