#!/usr/bin/env python3
# Simple test of the "equation trigger" functionality


import matplotlib.pyplot as plt
import numpy as np
from DREAM import DREAMSettings, runiface

import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


nt = 100
tMax = 1e-2

ds = DREAMSettings()

# Temperature equation
ds.eqsys.T_cold.setPrescribedData(100)

# Trigger alternative equation after 10% of the simulation
ds.eqsys.T_cold.trigger.setTimeTrigger(0.1 * tMax)
# Set self-consistent temperature
ds.eqsys.T_cold.trigger.equation.setType(T_cold.TYPE_SELFCONSISTENT)
ds.eqsys.T_cold.trigger.equation.setInitialProfile(100)

# Radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setWallRadius(0.30)
ds.radialgrid.setNr(20)

ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Add ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=1e19)

# Electric field
ds.eqsys.E_field.setPrescribedData(0.1)

# Set time step
ds.timestep.setNt(nt)
ds.timestep.setTmax(tMax)

ds.solver.setType(Solver.NONLINEAR)

ds.save('settings.h5')
do = runiface(ds, 'output.h5')

do.eqsys.T_cold.plot(r=range(0,20,2))
plt.show()

