#!/usr/bin/env python3
import h5py
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../py')


from DREAM.DREAMOutput import DREAMOutput

# Load simulation results
do = DREAMOutput("output_nbi.h5")

times = do.grid.t
r = do.grid.r



plt.scatter(do.grid.r, do.other.fluid.Tcold_NBI[-1, :])  # plot at final time
plt.ylabel('NBI deposition profile in DREAM')
plt.xlabel('r [cm]')
plt.show()





