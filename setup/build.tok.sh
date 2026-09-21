#!/bin/bash
#
# This script builds DREAM on Max-Planck IPP TOK systems.
source environment.tok.sh

cd $DREAMPATH && rm -rf build && mkdir build && cd build &&
cmake .. -DDREAM_BUILD_PYFACE=OFF -DPETSC_EXECUTABLE_RUNS=YES -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH
make -j8


