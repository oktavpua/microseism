#!/bin/bash
set -e

src=$1
if [ -z ${src} ]; then
    echo "Usage: $0 <microseism src>"
    exit 1
fi

if [ ! -d vexcl ]
then
	git clone https://github.com/ddemidov/vexcl
fi

CC=mpicc CXX=mpic++ cmake $src -DCMAKE_BUILD_TYPE=Release -DVEXCL_ROOT=${PWD}/vexcl -B${src}/build
cmake --build ${src}/build -- -j $(nproc)
