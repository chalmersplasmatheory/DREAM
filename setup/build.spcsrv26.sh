#!/bin/bash
#
# This script builds DREAM on SPCSRV26 (situated at the Swiss Plasma Center,
# Ecole Polytechnique Federale de Lausanne, Switzerland).
# The following environment variables can be specified (defaults indicated):
#
#   DREAMPATH=~/DREAM            (directory where DREAM will be looked for)
#
# If these variables are not already set in the environment, the
# "environment.spcsrv26.sh" script will initialize them to their default values.
#
# NOTE: When running DREAM, you will want to source the "environment.spcsrv26.sh"
# script.
#
# -------
# Written by Mathias Hoppe (2022-03-24)
# 
# CHANGELOG
#
# 2022-03-24
# - Created script.
#

source environment.spcsrv26.sh

function install_dream {
	echo $DREAMPATH
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	#cmake .. -DPETSC_EXECUTABLE_RUNS=YES &&
	cmake .. &&
	make -j64
}

install_dream

