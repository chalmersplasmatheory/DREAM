#!/usr/bin/bash
#
# This script builds DREAM on the Mahti system (CSC, Finland).
# environment variables can be specified (defaults indicated):
#
#   PETSC_DIR=~/petsc            (directory where PETSc will be looked for/downloaded to)
#   PETSC_ARCH=arch-linux-c-opt  (name of PETSc build to use)
#   DREAMPATH=~/DREAM            (directory where DREAM will be looked for)
#
# -------
# Written by Mathias Hoppe (2021-04-19)
#

source environment.mahti.sh

# Build PETSc
function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"
	# Configure with Intel MKL
	./configure --with-debugging=0 --COPTFLAGS="-O2 -march=native" --CXXOPTFLAGS="-O2 -march=native" -FOPTFLAGS="-O2 -march=native" CC=mpicc CXX=mpic++ FC=mpifort &&
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}

function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES &&
	make -j8
}

# Only build PETSc if this hasn't been done before
if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
    install_petsc
fi

install_dream

