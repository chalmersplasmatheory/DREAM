#!/usr/bin/bash
#
# Environment setup for the Puhti system (CSC, Finland)
#

# Set PETSc path (can be overridden by settings these variables explicitly)
if [ -z "$PETSC_DIR" ]; then
	export PETSC_DIR=~/petsc
fi
if [ -z "$PETSC_ARCH" ]; then
	export PETSC_ARCH=arch-linux-c-opt
fi
if [ -z "$DREAMPATH" ]; then
    export DREAMPATH=~/DREAM
fi

module purge
module load gcc openmpi
module load cmake gsl hdf5 openblas
module load python-env

#export HDF5_ROOT=$HDF5_INSTALL_ROOT


