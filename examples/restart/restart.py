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
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions


##############################
# PART I
##############################
ds1 = DREAMSettings()

E = 0.3     # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 1e3     # Temperature (eV)

# Set E_field
ds1.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds1.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds1.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Hot-tail grid settings
pmax = 2
ds1.hottailgrid.setNxi(10)
ds1.hottailgrid.setNp(500)
ds1.hottailgrid.setPmax(pmax)

ds1.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

# Set initial hot electron Maxwellian
ds1.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)

# Disable runaway grid
ds1.runawaygrid.setEnabled(False)

# Set up radial grid
ds1.radialgrid.setB0(5)
ds1.radialgrid.setMinorRadius(0.22)
ds1.radialgrid.setWallRadius(0.22)
ds1.radialgrid.setNr(1)

# Use the linear solver
ds1.solver.setType(Solver.LINEAR_IMPLICIT)

ds1.other.include('fluid/runawayRate')

# Set time stepper
ds1.timestep.setTmax(2e-1)
ds1.timestep.setNt(10)

# Set name of output file
ds1.output.setFilename('output1.h5')

# Save settings to HDF5 file
ds1.save('dream_settings_1.h5')


##############################
# PART II
##############################
ds2 = DREAMSettings(ds1)

ds2.fromOutput('output1.h5', ignore=['n_i'])
ds2.output.setFilename('output2.h5')
ds2.save('dream_settings_2.h5')

