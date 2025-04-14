#!/usr/bin/env python3
import h5py
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../py')

from DREAM.DREAMOutput import DREAMOutput

# Load simulation results
do = DREAMOutput("output_nbi.h5")

# --- Plot T_cold as function of radius at final time ---
times = do.grid.t
r = do.grid.r

T_cold = do.eqsys.T_cold[:,:]  # shape: [time, radius]
print("T_cold shape:", T_cold.shape)
print("T_cold final timestep:", T_cold[-1, :])

plt.figure()
plt.plot(r/1e-2, T_cold[-1, :], label='T_cold (final)')
plt.xlabel('r [cm]')
plt.ylabel('T_cold [eV]')
plt.title('Cold electron temperature at final time')
plt.grid(True)
plt.legend()
plt.show()

# --- Plot T_cold at center as function of time ---
plt.figure()
plt.plot(times, T_cold[:, 0], label='T_cold (r=0)')
plt.xlabel('Time [s]')
plt.ylabel('T_cold [eV]')
plt.title('Central T_cold vs Time')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

plt.scatter(do.grid.r, do.other.fluid.Tcold_NBI[-1, :])  # plot at final time
plt.ylabel('NBI')
print(do.grid.r)

plt.xlabel('r [cm]')
plt.show()
