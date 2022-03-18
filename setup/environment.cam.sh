#!/bin/bash
#
# Set PETSc path (can be overridden by settings these variables explicitly)
if [ -z "$PETSC_DIR" ]; then
	export PETSC_DIR=~/rds/rds-ukaea-ap001/RE/petsc/
fi
if [ -z "$PETSC_ARCH" ]; then
	export PETSC_ARCH=arch-linux-c-opt
fi
# Set DREAM path (based on location of this script)
if [ -z "$DREAMPATH" ]; then
    export DREAMPATH=~/rds/rds-ukaea-ap001/RE/DREAM
fi

# Modules needed for building and running DREAM
# Clean up module environment
module purge

# Load required modules
module load gcc/11 slurm
module load hdf5-1.10.1-gcc-5.4.0-z7bre2k
module load cmake gsl/2.4
module load python

#Active Python conda environment (needed for h5py)
module load miniconda/3
cd ~/rds/rds-ukaea-ap001/RE/
conda init
conda deactivate
conda activate DREAM_PYTHON
cd -
