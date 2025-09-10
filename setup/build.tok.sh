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
# 2024-03-25
# - Minor fix, it can now be installed on TOK machines. Does not seem to work with MKL, however... /Peter
#
source environment.tok.sh

function install_petsc {

	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"

	# Configure with Intel MKL
	./configure --with-debugging=0 --COPTFLAGS="-O3 -march=native -mtune=native" 
--CXXOPTFLAGS="-O3 -march=native -mtune=native" --FOPTFLAGS="-O3 --march=native -mtune=native" 
--with-mkl_pardiso=1 --with-mkl_pardiso-dir=$MKL_HOME --with-blaslapack-dir=$MKL_HOME -with-mpi=0

	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}
function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DDREAM_BUILD_PYFACE=OFF -DPETSC_EXECUTABLE_RUNS=YES -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH
	make -j8
}

install_petsc
install_dream


