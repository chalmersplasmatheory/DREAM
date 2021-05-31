#!/bin/bash
#
# This script builds PETSc and DREAM on Max-Planck IPP TOK systems.
# The following environment variables can be specified (defaults indicated):
#
#   DREAMPATH=~/DREAM            (directory where DREAM will be looked for)
#
# If these variables are not already set in the environment, the
# "environment.tok.sh" script will initialize them to their default values.
#
# NOTE: When running DREAM, you will want to source the "environment.tok.sh"
# script.
#
# -------
# Written by Mathias Hoppe (2021-04-01)
#

source environment.tok.sh

function install_dream {
	cd "$DREAM_DIR" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES &&
	make -j8
}

install_dream

