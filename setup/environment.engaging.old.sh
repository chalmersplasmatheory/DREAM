#!/bin/sh
#
# Environment setup for the engaging cluster (PSFC/MIT)

# Set DREAM path (based on location of this script)
if [ -z "$DREAMPATH" ]; then
	export DREAMPATH=$(realpath "$(dirname $0)"/..)
fi

# Clean up module environment
module purge

# Load required modules
module load cmake/3.17.3
#module load gcc/11.2.0
module load intel/2021.3.0

module load gsl/2.5
export LIB=/home/software/gsl/2.5/lib:${LD_LIBRARY_PATH}
export INCLUDE=/home/software/gsl/2.5/include:${CPATH}

module load hdf5/1.10.5-cxx

module load python/3.9.4

