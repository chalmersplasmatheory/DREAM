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
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

ds = DREAMSettings()

times  = [0]
radius = [0, 1]

E_selfconsistent = 1
T_selfconsistent = 0

# Set E_field 
if not E_selfconsistent:
    efield = 2000*np.ones((len(times), len(radius)))
    ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
else:
    ds.eqsys.E_field = ElectricField(Efield.TYPE_SELFCONSISTENT, efield=0.0)
    ds.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = 5000)

if not T_selfconsistent:
    temperature = 10 * np.ones((len(times), len(radius)))
    ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)
else:
    ds.eqsys.T_cold = ColdElectronTemperature(ttype=T_cold.TYPE_SELFCONSISTENT, temperature=10.0)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=1e20)

# Hot-tail grid settings
#pmax = 0.1
#ds.hottailgrid.setNxi(30)
#ds.hottailgrid.setNp(500)
pmax = 0.05
ds.hottailgrid.setNxi(4)
ds.hottailgrid.setNp(100)
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
ds.radialgrid.setMinorRadius(1)
ds.radialgrid.setNr(10)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)

#ds.other.include('nu_s')
#ds.other.include('all')
ds.other.include('fluid')


# Set time stepper
#ds.timestep.setTmax(1e-2)
#ds.timestep.setNt(100)
ds.timestep.setTmax(.0001)
ds.timestep.setNt(10)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

