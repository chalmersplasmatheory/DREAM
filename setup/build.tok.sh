#!/bin/bash
#
# This script builds PETSc and DREAM on Max-Planck IPP TOK systems.
# The following environment variables can be specified (defaults indicated):
#
#   DREAMPATH=~/DREAM            (directory where DREAM will be looked for)
#
# If these variables are not already set in the environment, the
# "environment.tok.sh" script will initialize them to their default values.
#
# NOTE: When running DREAM, you will want to source the "environment.tok.sh"
# script.
#
# -------
# Written by Mathias Hoppe (2021-04-01)
#
# CHANGELOG
#
# 2021-08-09
# - Updated to build a custom version of PETSc with Intel MKL
#

source environment.tok.sh

function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"
	# Configure with Intel MKL
    ./configure --with-debugging=0 --COPTFLAGS="-O3 -march=native -mtune=native" --CXXOPTFLAGS="-O3 -march=native -mtune=native" --FOPTFLAGS="-O3 -march=native -mtune=native" --with-mkl_pardiso=1 &&
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}
function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES &&
	make -j8
}

install_petsc
install_dream
