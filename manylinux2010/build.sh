#!/bin/bash -ve

HERE=`dirname "$0"`
cd $HERE/..

for i in 36 37 38 39 310; do
    docker build -f manylinux2010/wheel${i}.docker . -t pybdsf${i}
    docker run -v `pwd`/manylinux2010:/manylinux2010 pybdsf${i} sh -c "cp /output/*.whl /manylinux2010/."
done

