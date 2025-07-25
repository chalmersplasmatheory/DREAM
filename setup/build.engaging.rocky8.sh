#!/bin/sh
#
# Build DREAM on engaging (rocky8 machine) (PSFC/MIT)

#source environment.engaging.rocky8.sh
THISPATH=$(realpath "$(dirname $0)")
source "$THISPATH/environment.engaging.rocky8.sh"

if [ -z "$SLURM_JOB_ID" ]; then
	echo "WARNING: No SLURM job allocation has been detected. Consider running"
	echo "WARNING:   salloc -p mit_normal -c 16"
	echo "WARNING: before proceeding with the installation to significantly improve performance."
	echo "WARNING: (now waiting for 5 seconds...)"
	sleep 5s
fi

function install_petsc {
	if [ ! -d "$PETSC_DIR" ]; then
		git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
	fi

	cd "$PETSC_DIR"

	./configure PETSC_ARCH=$PETSC_ARCH --with-mpi=0 --download-make \
		--with-debugging=0 --COPTFLAGS="-O3 -march=native -mtune=native" \
	       	--CXXOPTFLAGS="-O3 -march=native -mtune=native" \
		--FOPTFLAGS="-O3 -march=native -mtune=native" \
		--with-mkl_pardiso=1 \
		--with-mkl_pardiso-dir=/orcd/software/community/001/rocky8/intel/2024.2.1/mkl/2024.2 \
		--with-blaslapack-dir=/orcd/software/community/001/rocky8/intel/2024.2.1/mkl/2024.2 \
		--with-cc=icx --with-cxx=icpx --with-fc=ifx
	
	make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH
}

function install_hdf5 {
	if [ ! -d "$HDF5_BUILD_DIR" ]; then
		mkdir -p "$HDF5_BUILD_DIR"
		cd "$HDF5_BUILD_DIR"
		wget https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
		tar -xzf hdf5-1.14.6.tar.gz
		cd hdf5-1.14.6
	else
		cd "$HDF5_BUILD_DIR/hdf5-1.14.6"
	fi

	./configure --prefix="$HDF5_DIR" --enable-cxx --enable-hl CC=icx CXX=icpx
	make
	make install
}

function install_dream {
	cd "$DREAMPATH" && rm -rf build && mkdir build && cd build &&
	cmake .. -DPETSC_EXECUTABLE_RUNS=YES \
		-DCMAKE_CXX_COMPILER=$HDF5_DIR/bin/h5c++ \
		-DGSL_ROOT_DIR=/home/software/modulefiles/gsl/2.5
	
	make -j 8
}

function is_python_installed {
	python -c "import h5py, matplotlib, numpy, packaging, scipy"
	return $?
}

if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
	install_petsc
fi

if [ ! -d "$HDF5_DIR" ]; then
	install_hdf5
fi

is_python_installed
if [ $? -ne 0 ]; then
	echo "ERROR: Unable to find required Python packages. Are you really running with the appropriate Conda environment?"
	exit -1
fi

install_dream
