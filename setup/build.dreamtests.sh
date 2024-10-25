#!/bin/bash
#
# This script builds DREAM on LOTTA (situated at KTH, Sweden).
# The following environment variables can be specified (defaults indicated):
#
#   DREAMPATH=~/DREAM            (directory where DREAM will be looked for)
#
# If these variables are not already set in the environment, the
# "environment.lotta.sh" script will initialize them to their default values.
#
# NOTE: When running DREAM, you will want to source the "environment.lotta.sh"
# script.
#
# -------
# Written by Mathias Hoppe (2023-04-04)
# 
# CHANGELOG
#
# 2023-04-04
# - Created script.
#

source environment.lotta.sh

# Build PETSc
function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"

	# Configure with Intel MKL?
	./configure --with-debugging=0 --COPTFLAGS=-O2 --CXXOPTFLAGS=-O2 --FOPTFLAGS=-O2 --with-mkl_pardiso=1 --with-mkl_pardiso-dir=/opt/intel/oneapi/mkl/latest --with-blaslapack-dir=/opt/intel/oneapi/mkl/latest --with-mpi=0 &&
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}

function install_dream {
	echo $DREAMPATH
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH -DPETSC_EXECUTABLE_RUNS=YES &&
	make -j8
}

HAS_PETSC=1
if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
	install_petsc

	if [ $? -ne 0 ]; then
		echo "ERROR: Failed to install PETSc"
		HAS_PETSC=0
	fi
fi

if [ $HAS_PETSC -gt 0 ]; then
	install_dream
fi

