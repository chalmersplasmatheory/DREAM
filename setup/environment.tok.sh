#!/bin/sh
#
# Modules needed for building and running DREAM

# Set PETSc path (can be overridden by settings these variables explicitly)
if [ -z "$PETSC_DIR" ]; then
	export PETSC_DIR=~/petsc
fi
if [ -z "$PETSC_ARCH" ]; then
	export PETSC_ARCH=arch-linux-c-opt
fi
# Set DREAM path (based on location of this script)
if [ -z "$DREAMPATH" ]; then
    export DREAMPATH=$(realpath "$(dirname $0)"/..)
fi

# Clean up module environment
module purge

# Load required modules
module load gcc/10
module load anaconda/3/2020.02
module load hdf5-serial/1.12.0 openmpi
module load mkl/2020.4
module load cmake git gsl

alias dreamviz="python -i $DREAMPATH/py/cli/cli.py"

