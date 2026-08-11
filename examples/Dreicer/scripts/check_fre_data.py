#!/usr/bin/env python3
#
# Check if f_re data exists and examine its structure
#

import numpy as np
import sys

sys.path.append('../../py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput

# Load data
do = DREAMOutput('../outputs/dreicer_with_fre_output.h5')

print("Available quantities in eqsys:")
print(dir(do.eqsys))

print("\nChecking for f_re:")
if hasattr(do.eqsys, 'f_re'):
    print("f_re exists")
    print(f"f_re shape: {do.eqsys.f_re[:].shape}")
else:
    print("f_re does not exist")

print("\nChecking for f_hot:")
if hasattr(do.eqsys, 'f_hot'):
    print("f_hot exists")
    print(f"f_hot shape: {do.eqsys.f_hot[:].shape}")

print("\nRunaway electron density:")
if hasattr(do.eqsys, 'n_re'):
    n_re = do.eqsys.n_re[:]
    print(f"n_re shape: {n_re.shape}")
    print(f"n_re values: {n_re.flatten()}")

# Close the file
do.close()