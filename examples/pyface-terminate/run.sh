#!/usr/bin/env bash

# Preload Intel MKL
export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_def.so:/opt/intel/mkl/lib/intel64/libmkl_avx2.so:/opt/intel/mkl/lib/intel64/libmkl_core.so:/opt/intel/mkl/lib/intel64/libmkl_intel_lp64.so:/opt/intel/mkl/lib/intel64/libmkl_intel_thread.so:/opt/intel/mkl/lib/intel64_lin/libiomp5.so

# Run DREAM through pyface
if [[ $# -eq 0 ]]; then
    ./generate.py
elif [[ $# -eq 1 ]]; then
    python3 "$1"
else
    echo "ERROR: Too many command line arguments provided"
    exit 1
fi

