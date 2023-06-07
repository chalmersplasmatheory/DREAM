#!/usr/bin/env python3
#
# Set up a simulation to test the ion source.
#

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM import DREAMSettings, runiface
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


ds = DREAMSettings()

E  = 0.1
n0 = 5e19
T = 800

tMax = 0.1
nt = 100

# Configure physics
ds.eqsys.E_field.setPrescribedData(E)
ds.eqsys.T_cold.setPrescribedData(T)

# Configure ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, Z0=1, n=n0)
ds.eqsys.n_i.addIonSource('D', dNdt=1e19)  # constant Z0=0 source

ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.42)
ds.radialgrid.setWallRadius(0.44)
ds.radialgrid.setNr(10)

ds.solver.setType(Solver.LINEAR_IMPLICIT)

ds.other.include('fluid')

ds.timestep.setTmax(tMax)
ds.timestep.setNt(nt)

ds.output.setTiming(stdout=True, file=True)

ds.save('settings.h5')

do = runiface(ds, 'output.h5')

#####################
# PLOT RESULTS
#####################
fig, axs = plt.subplots(1, 2, figsize=(12, 3))

do.eqsys.n_i['D'][1].plot(t=[0, 20, 40, 60, 80, 100], ax=axs[0])

n1 = do.eqsys.n_i['D'][1]
N  = n1.integral()
dn = np.diff(n1[:,-1]) / np.diff(do.grid.t)
dN = np.diff(N[:]) / np.diff(do.grid.t)

#axs[1].plot(do.grid.t[1:], dn)
#axs[1].set_ylim([0, 1.1*max(dn)])
axs[1].plot(do.grid.t[1:], dN)
axs[1].set_ylim([0, 1.1*max(dN)])
axs[1].set_xlim([0, do.grid.t[-1]])
axs[1].set_title('Source magnitude')

print(f'dN/dt = {dN[-1]} particles/s')

fig.tight_layout()
plt.show()

