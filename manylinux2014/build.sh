#!/bin/bash
#
# Script to make python wheels for several versions

set -eu

SCRIPT_DIR=$(cd $(dirname $0) && pwd)
ROOT_DIR=$(git rev-parse --show-toplevel)

for i in 36 37 38 39 310; do
    docker build \
        --file ${SCRIPT_DIR}/wheel${i}.docker \
        --tag pybdsf${i} \
        ${ROOT_DIR}
    dockerid=$(docker create pybdsf${i})
    docker cp ${dockerid}:/dist ${ROOT_DIR}
    docker rm ${dockerid}
done

