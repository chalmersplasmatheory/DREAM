
#!/bin/bash
#
# This script builds PETSc and DREAM on Cambridge CSD3 cluster (./build.cam.sh)
# (https://docs.hpc.cam.ac.uk/hpc/index.html)
#
# Please do ". environment.cam.sh" before running ./build.cam.sh (otherwise the initialization
# of conda fails).
# The "environment.cam.sh" script will initialize DREAMPATH to its default value, and load
# the required modules. petsc is installed locally in the shared RE directory, if you want
# to install your own version of petsc and DREAM, you can modify the paths in environment.cam.sh.
# If left empty, it will download and compile it in your home directory.
#
# NOTE: When running DREAM, you will also want to source the "environment.cam.sh"
# script beforehand.
#
# -------
# Written by Alexandre Fil (2022-03-01)
# 

function install_petsc {
	if [ -d "$PETSC_DIR" ]
	then
	    if [ "$(ls -A $PETSC_DIR)" ]; then
		echo "PETSC is already installed in $PETSC_DIR"
	    else
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
		cd "$PETSC_DIR"
		./configure --download-f2cblaslapack=1 --with-debugging=0 --COPTFLAGS="-O3 -march=native -mtune=native" --CXXOPTFLAGS="-O3 -march=native -mtune=native" --FOPTFLAGS="-O3 -march=native -mtune=native" --download-superlu --download-mpich=1 &&
		make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
	    fi
	fi
}

function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES -DGSL_ROOT_DIR=/usr/local/Cluster-Apps/gsl/2.4/ -DMPI_CXX_COMPILER=~/rds/rds-ukaea-ap001/RE/petsc/arch-linux-c-opt/bin/ -DPETSC_DIR="$PETSC_DIR" -DPETSC_ARCH="$PETSC_ARCH" &&
	make -j8
}

install_petsc
install_dream

