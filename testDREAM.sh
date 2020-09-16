#!/usr/bin/env bash
#
# This script run all tests available for DREAM and reports whether
# the tests pass or not.
#
# THIS SCRIPT SHOULD BE RUN BEFORE SUBMITTING A PULL REQUEST.
# --------------------------------------------------------------
# RUN THIS SCRIPT FROM THE DREAM ROOT DIRECTORY:
# 
#   $ ./testDREAM.sh
#
# ##############

function logerror {
    >&2 echo -e "\x1B[1;31mERROR \x1B[33m$1\x1B[31m: $2\x1B[0m"
    exit $1
}

# Make sure the kernel has been rebuilt
NCPUS=$(grep -c ^processor /proc/cpuinfo)
cd build && make -j "$NCPUS" && cd ..

if [ $? -ne 0 ]; then
    logerror 1 "Failed to build DREAM."
fi

# Run C++ tests
build/tests/cxx/dreamtests all

if [ $? -ne 0 ]; then
    logerror 2 "C++ tests failed."
fi

# Run physics/Python tests
tests/physics/runtests.py all

if [ $? -ne 0 ]; then
    logerror 3 "Physics/Python tests failed."
fi

