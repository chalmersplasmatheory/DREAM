#!/bin/sh
#
# Modules needed for building and running DREAM

export HDF5_USE_FILE_LOCKING=FALSE

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

# Clean up module environment
module purge

# Load required modules
module load PrgEnv-intel
#module swap gcc gcc/10.1.0
module load cmake/3.18.2 gsl
module load python/3.8-anaconda-2020.11

# PETSc dependencies
module load cray-hdf5
module load impi

# Set up Python environment (see "build.cori.sh")
source activate dreamenv

alias dreamviz="python -i $DREAMPATH/py/cli/cli.py"

