#!/bin/sh
#
# Modules needed for building and running DREAM
#
# CHANGELOG
#
# 2024-03-25
# - Updated modules. /Peter

# Set PETSc path (can be overridden by settings these variables explicitly)
if [ -z "$PETSC_DIR" ]; then
	export PETSC_DIR=~/petsc
fi
if [ -z "$PETSC_ARCH" ]; then
	export PETSC_ARCH=arch-linux-c-opt
fi
# Set DREAM path (based on location of this script)
if [ -z "$DREAMPATH" ]; then
    export DREAMPATH="$(dirname -- "$(realpath -- "$0")")"
	echo $DREAMPATH
fi

# Clean up module environment
module purge

# Load required modules
module load gcc/13
module load anaconda/3/2023.03 
module load hdf5-serial/1.14.1 openmpi/4.1
module load mkl/2023.1
module load cmake git gsl/2.4

alias dreamviz="python -i $DREAMPATH/py/cli/cli.py"

