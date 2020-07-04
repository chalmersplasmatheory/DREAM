#!/usr/bin/env python3
#
# This example shows how to set up a simple CODE-like runaway
# scenario in DREAM. The simulation uses a constant temperature,
# density and electric field, and generates a runaway current
# through the electric field acceleration.
#
# Run as
#
#   $ ./basic.py
#   $ ../../build/iface/dreami dream_settings.h5
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions


ds = DREAMSettings()

#E = 0.3     # Electric field strength (V/m)
E = 6.745459970079014
n = 5e19    # Electron density (m^-3)
#T = 1e3     # Temperature (eV)
T = 100

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.avalanche = False

# Hot-tail grid settings
pmax = 2
ds.hottailgrid.setNxi(10)
ds.hottailgrid.setNp(500)
ds.hottailgrid.setPmax(pmax)

ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(1)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.solver.setTiming(True)

ds.other.include('fluid/runawayRate')

# Set time stepper
ds.timestep.setTmax(2e-1)
ds.timestep.setNt(10)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

