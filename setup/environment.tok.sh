#!/bin/sh
#
# Modules needed for building and running DREAM

# Set PETSc path (can be overridden by settings these variables explicitly)
if [ -z "$DREAMPATH" ]; then
    export DREAMPATH=~/DREAM
fi

# Clean up module environment
module purge

# Load required modules
module load gcc/10
module load anaconda/3/2020.02
module load hdf5-serial/1.12.0 openmpi
module load petsc-real
module load cmake git gsl

alias dreamviz="python -i $DREAMPATH/py/cli/cli.py"

