#!/bin/bash

# This is a script to build the modules and run the test suite in the base
# Docker container.

die() {
  echo "Error: $@" 1>&2
  exit 1;
}

cd /usr/src/ITKMinimalPathExtraction-build || die "Could not cd into the build directory"

cmake \
  -G Ninja \
  -DITK_DIR:PATH=/usr/src/ITK-build \
  -DCMAKE_BUILD_TYPE:STRING=Release \
    /usr/src/ITKMinimalPathExtraction || die "CMake configuration failed"
ctest -VV -D Experimental || die "ctest failed"

examples_build_dir=/usr/src/ITKMinimalPathExtraction-build/examples
mkdir -p $examples_build_dir
cd $examples_build_dir
cmake \
  -G Ninja \
  -DITK_DIR:PATH=/usr/src/ITK-build \
  -DCMAKE_BUILD_TYPE:STRING=Release \
    /usr/src/ITKMinimalPathExtraction/examples || die "Example CMake configuration failed"
ninja || die "examples build failed"
