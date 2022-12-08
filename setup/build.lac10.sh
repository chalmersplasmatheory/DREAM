#!/usr/bin/bash
#
# This script builds DREAM on LAC10 (situated at the Swiss Plasma Center,
# Ecole Polytechnique Federale de Lausanne, Switzerland). Note that the DREAM is
# only able to run on LAC10, since the required compilers are not avaiable on
# other LAC machines (as of 2022.12.08).
#
# The following environment variables can be specified (defaults indicated):
#
#   DREAMPATH=~/DREAM              # (directory where DREAM will be installed)
#
# If these variables are not already set in the environment, the
# "environment.lac10.sh" script will initialize them to their default values.
#
# NOTE: When running DREAM, you will want to source the "environment.lac10.sh"
# script.
#
# -------
# Written by Mathias Hoppe (2022-12-08)
#
# CHANGELOG
#
# 2022-12-08
# - Created script.
#

source environment.lac10.sh

# Build PETSc
function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"

	# Configure with Intel MKL?
	./configure --with-debugging=0 --COPTFLAGS=-O2 --CXXOPTFLAGS=-O2 --FOPTFLAGS=-O2 --with-blas-lib=/usr/lib64/libblas.so.3.8.0 --with-lapack-lib=/usr/lib64/liblapack.so.3.8.0 --with-mpi=0 &&
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}

# Install Python packages
function is_python_installed {
	python3 -c "import h5py, matplotlib, numpy, packaging, scipy"
	return $?
}
function install_python {
	echo "Installing required Python packages"
	pip install --user h5py matplotlib numpy packaging scipy
}

function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH &&
	make -j8
}

if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
	install_petsc
fi

is_python_installed
if [ $? -ne 0 ]; then
	install_python
fi

install_dream

