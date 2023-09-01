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

tMax = 1e-4
nt = 100

# Configure physics
ds.eqsys.E_field.setPrescribedData(E)
ds.eqsys.T_cold.setPrescribedData(T)

# Configure ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, Z0=1, n=n0)
ds.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC, n=1e19, init_equil=True)

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
fig, axs = plt.subplots(1, 1, figsize=(6, 3))

do.eqsys.n_i['Ne'].plot(ax=axs)

fig.tight_layout()
plt.show()

