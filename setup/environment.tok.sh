#!/bin/sh
#
# Modules needed for building and running DREAM
#
# CHANGELOG
#
# 2024-03-25
# - Updated modules. /Peter
#
# 2025-07-02
# - Updated modules again. /Peter
#
# 2025-02-06
# - Difficulties with setting up PETSc. Now uses available module installation. 
#	Note that the module does not use the PETSC_ARCH split directory layout.
# 	/Peter
#	
# 2026-09-21
# - Updated modules. /Peter (and Victor Svensson)

if [[ -z "${DREAMPATH:-}" ]]; then
	export DREAMPATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)" 
fi



# Clean up module environment
module purge

# Load required modules
module load gcc
module load hdf5-serial/1.14.1 openmpi/5.0
module load petsc-real-double/3.25	# load PETSc module rather than install
module load cmake git gsl
module load python-waterboa

alias dreamviz="python -i $DREAMPATH/py/cli/cli.py"
alias dreami="$DREAMPATH/build/iface/dreami"
