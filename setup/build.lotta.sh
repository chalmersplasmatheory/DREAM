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
	#./configure --with-debugging=0 --COPTFLAGS=-O2 --CXXOPTFLAGS=-O2 --FOPTFLAGS=-O2 --with-blas-lib=/usr/lib/x86_64-linux-gnu/blas64/libblas64.so.3.10.0 --with-lapack-lib=/usr/lib/x86_64-linux-gnu/lapack64/liblapack64.so.3.10.0 --with-mpi=0 &&
	./configure --with-debugging=0 --COPTFLAGS=-O2 --CXXOPTFLAGS=-O2 --FOPTFLAGS=-O2 --with-blas-lib=/usr/lib/x86_64-linux-gnu/blas/libblas.so --with-lapack-lib=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so  --with-mpi=0 &&
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}

function install_dream {
	echo $DREAMPATH
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH &&
	make -j64
}

HAS_PETSC=1
if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
	install_petsc

	if [ -$? -ne 0 ]; then
		echo "ERROR: Failed to install PETSc"
		HAS_PETSC=0
		#rm -rf "$PETSC_DIR/$PETSC_ARCH"
	fi
fi

if [ $HAS_PETSC -gt 0 ]; then
	install_dream
fi

