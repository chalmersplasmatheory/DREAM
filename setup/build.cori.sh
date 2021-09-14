#!/bin/bash
#
# This script builds PETSc and DREAM on the NERSC Cori system. The following
# environment variables can be specified (defaults indicated):
#
#   PETSC_DIR=~/petsc            (directory where PETSc will be looked for/downloaded to)
#   PETSC_ARCH=arch-linux-c-opt  (name of PETSc build to use)
#   DREAM_DIR=~/DREAM            (directory where DREAM will be looked for)
#
# If these variables are not already set in the environment, the
# "environment.cori.sh" script will initialize them to their default values.
#
# NOTE: When running DREAM, you will want to source the "environment.cori.sh"
# script.
#
# A number of different compilers are available on Cori, but several
# combinations of software have been observed to fail for various reasons. This
# configuration uses Intel compilers, Intel MPI and a custom build of PETSc. It
# is probably possible to instead use Cray compilers and MPICH, although this
# has not been tried. GNU compilers fail to compile DREAM because of a
# configuration error with the HDF5 library.
#
# PYTHON
# In addition to loading the Python module, it is necessary to set up a custom
# Python environment. This is done once by running the following after loading
# the python module:
#
#   conda create --name dreamenv python=3.8
#
# Then you should run
#
#   source activate dreamenv
#
# once in every shell where you would like to use Python (e.g. to plot output
# using the DREAM command-line interface). This command is also included in the
# "environment.cori.sh" script. Note also that the required Python packages must
# be installed manually after creating the Python environment above. At the time
# of writing this, the full command required should be
#
#   conda install h5py matplotlib numpy packaging scipy
#   
# -------
# Written by Mathias Hoppe (2021-04-01)
#

source environment.cori.sh

# Build PETSc
function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"
	# Configure with Intel MKL
	./configure --with-debugging=0 --COPTFLAGS=-O2 --CXXOPTFLAGS=-O2 -FOPTFLAGS=-O2 --with-mkl_pardiso=1 CC=mpiicc CXX=mpiicpc FC=mpiifort --with-mkl_pardiso-dir=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl --with-blaslapack-dir=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl &&
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
}

function install_dream {
	cd "$DREAM_DIR" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH &&
	make -j8
}

# Only build PETSc if this hasn't been done before
if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
    install_petsc
fi

install_dream

