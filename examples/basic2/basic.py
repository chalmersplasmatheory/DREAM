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
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ElectricField as Efield

from DREAM.Settings.Equations.ElectricField import ElectricField

ds = DREAMSettings()

times  = [0]
radius = [0, 1]

# Set E_field 
efield = 2000*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)

# Set self-consistent E-field evolution

#ds.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
#ds.eqsys.E_field = ElectricField(Efield.TYPE_SELFCONSISTENT, efield=1.0)

# Set n_cold (prescribed; it is automatically calculated self-consistently otherwise)
#density = 1e20 * np.ones((len(times), len(radius)))
#ds.eqsys.n_cold.setPrescribedData(density=density, times=times, radius=radius)

# Set temperature
temperature = 10 * np.ones((len(times), len(radius)))
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=1e20)

# Hot-tail grid settings
#pmax = 0.1
#ds.hottailgrid.setNxi(30)
#ds.hottailgrid.setNp(500)
pmax = 0.1
ds.hottailgrid.setNxi(5)
ds.hottailgrid.setNp(600)
ds.hottailgrid.setPmax(pmax)


#ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_ULTRA_RELATIVISTIC
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_NON_SCREENED
#ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
#ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_NEGLECT
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT

# Set initial Maxwellian @ T = 1 keV, n = 5e19, uniform in radius
#ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=5e19, rT0=0, T0=1e3)
ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=1e20, rT0=0, T0=10)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setWallRadius(0.22)
ds.radialgrid.setNr(3)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Also output collision frequencies
# ('nu_s' stores ALL slowing-down frequencies; one can also specify
# each frequency separately:
#   hottail/nu_s, hottail/nu_s_fr, hottail/nu_s_f1, hottail/nu_s_f2,
#   runaway/nu_s, runaway/nu_s_fr, runaway/nu_s_f1, runaway/nu_s_f2
#ds.other.include('nu_s')
#ds.other.include('all')
ds.other.include('nu_s','nu_D','fluid')


# Set time stepper
#ds.timestep.setTmax(1e-2)
#ds.timestep.setNt(100)
ds.timestep.setTmax(1e-7)
ds.timestep.setNt(2)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

