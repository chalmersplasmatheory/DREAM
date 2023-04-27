#!/bin/sh
#
# Environment setup for LOTTA (@ KTH)

# Assume codes are installed directly under the user's home
USERHOME=$(eval echo ~$USER)

if [ -z "$PETSC_DIR" ]; then
	export PETSC_DIR="$USERHOME/petsc"
fi
if [ -z "$PETSC_ARCH" ]; then
	export PETSC_ARCH=arch-linux-c-opt
fi
if [ -z "$DREAMPATH" ]; then
	export DREAMPATH="$USERHOME/DREAM"
fi

alias dreamviz="python3 -i $DREAMPATH/py/cli/cli.py"

