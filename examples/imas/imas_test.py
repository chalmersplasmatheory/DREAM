#!/usr/bin/env python3
#
# A basic example script on the imas wrapper for DREAM. This script generates a basic
# DREAM input file which can be passed to 'dreami'. It is based on the basic.py example.
#
# Run as
#
#   $ ./imas_test.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.IMAS.IMAS as iw


ds = iw.readInIDSSlice(64614, 9999, 'tcv', log='test', setUpDream=True, wall_radius=0.25)

ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e19)

# Set fluid runaway models
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Disable hottail grid
ds.hottailgrid.setEnabled(False)

r = np.linspace(0, max(ds.radialgrid.r_f), len(ds.radialgrid.r_f))
rT = np.linspace(0, max(ds.radialgrid.r_f), len(ds.eqsys.T_cold.temperature))
n0 = 4e18*np.ones(len(ds.radialgrid.r_f))

# rn, n0, rT, T0 = ...  get profiles of density and temperature
ds.eqsys.f_hot.setInitialProfiles(rn0=r, n0=n0, rT0=rT, T0=ds.eqsys.T_cold.temperature)
ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_ANALYTIC_ALT_PC)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Use the linearly implicit solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Set time stepper
ds.timestep.setTmax(1.0e-5)
ds.timestep.setNt(400)

ds.other.include('fluid')

# Save settings to HDF5 file
ds.save('dream_settings.h5')
