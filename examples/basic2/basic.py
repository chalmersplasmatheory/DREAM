#!/usr/bin/env python3
#
# A very basic DREAM Python example. This script generates a basic
# DREAM input file which can be passed to 'dreami'.
#
# Run as
#
#   $ ./basic.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


ds = DREAMSettings()

times  = [0]
radius = [0, 1]

# Set E_field
efield = 1e-1*np.ones((len(times), len(radius)))
ds.equationsystem.E_field.setPrescribedData(efield=efield, times=times, radius=radius)

# Set n_cold (prescribed; it is automatically calculated self-consistently otherwise)
#density = 1e20 * np.ones((len(times), len(radius)))
#ds.equationsystem.n_cold.setPrescribedData(density=density, times=times, radius=radius)

# Set temperature
temperature = 1e4 * np.ones((len(times), len(radius)))
ds.equationsystem.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)

# Set ions
ds.equationsystem.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
ds.equationsystem.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=1e20)

# Hot-tail grid settings
pmax = 0.75
ds.hottailgrid.setNxi(10)
ds.hottailgrid.setNp(100)
ds.hottailgrid.setPmax(pmax)

# Set initial Maxwellian @ T = 1 keV, n = 5e19, uniform in radius
ds.equationsystem.f_hot.setInitialProfiles(rn0=0, n0=5e19, rT0=0, T0=1000)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(10)

# Use the linear solver
ds.solver.setType(Solver.TYPE_LINEAR_IMPLICIT)

# Set time stepper
#ds.timestep.setTmax(1e-2)
#ds.timestep.setNt(100)
ds.timestep.setTmax(1e-2)
ds.timestep.setNt(10)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

