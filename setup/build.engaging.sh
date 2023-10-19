#!/bin/bash
#
# Script for building DREAM on the "engaging" cluster (PSFC/MIT)
#

source environment.engaging.sh


function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"

	# Load MKL
	module load intel/2021.3.0

	#./configure PETSC_ARCH=$PETSC_ARCH --with-mpi=0 --download-fblaslapack=1 --download-make
	./configure PETSC_ARCH=$PETSC_ARCH --with-mpi=0 --download-make \
		--with-mkl_pardiso=1 \
		--with-mkl_pardiso-dir=/home/software/intel/2021.3.0/mkl/latest \
		--with-blaslapack-dir=/home/software/intel/2021.3.0/mkl/latest \
		CC=$CC CXX=$CXX F77=$F77
	make
}

function is_python_installed {
	python3 -c "import h5py, matplotlib, numpy, packaging, scipy"
	return $?
}

function install_python {
	pip install matplotlib packaging scipy
}

function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DCMAKE_INSTALL_PREFIX=${PWD} -DPETSC_EXECUTABLE_RUNS=YES -DPETSC_WITH_MPI=OFF -DGSL_ROOT_DIR=/home/software/modulefiles/gsl/2.5 &&
	make
}

if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
	install_petsc
fi

is_python_installed
if [ $? -ne 0 ]; then
	install_python
fi

install_dream

