#!/bin/sh
#
# Environment setup for the engaging cluster (rocky8 machine) (PSFC/MIT)

# Set DREAM path (based on location of this script)
if [ -z "$DREAMPATH" ]; then
	export DREAMPATH=$(realpath "$(dirname $0)"/..)
fi

# Clean up module environment
module purge

# Load required modules
module load community-modules intel/2024.2.1 gcc/12.2.0 cmake/3.27.9 gsl/2.7.1 miniforge/24.3.0-0

# Check if the 'dream' conda environment is available
dreamgrep=$(conda env list | grep dream)
if [ $? -ne 0 ]; then
	# Create new conda environment for DREAM
	conda create --name dream h5py matplotlib numpy packaging scipy
fi

conda activate dream

export PETSC_DIR=$HOME/packages/petsc
export PETSC_ARCH=arch-linux-c-opt
export HDF5_DIR=$HOME/packages/hdf5
export HDF5_BUILD_DIR=$HOME/packages/hdf5-build

