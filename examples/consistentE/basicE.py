#!/usr/bin/env python3
#
# This example shows how to set up a self-consistent fluid DREAM run,
# where no kinetic equations are solved, but the electric field and
# temperature are evolved self-consistently.
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

# set collision settings
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
#ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_NEGLECT
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
#ds.collisions.lnlambda = Collisions.LNLAMBDA_CONSTANT
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT

#############################
# Set simulation parameters #
#############################

B0 = 5          # magnetic field strength in Tesla
E_initial = 30  # initial electric field in V/m
E_wall = 20     # boundary electric field in V/m
T_initial = 10  # initial temperature in eV

Tmax = 1e-3     # simulation time in seconds
Nt = 3          # number of time steps
Nr = 5          # number of radial grid points
Np = 200        # number of momentum grid points
Nxi = 5         # number of pitch grid points
pMax = 0.04     # maximum momentum in m_e*c
times  = [0]    # times at which parameters are given
radius = [0, 1] # span of the radial grid
radius_wall = 2 # location of the wall 

hotTailGrid_enabled = False

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)
# Set time stepper
ds.timestep.setTmax(Tmax)
ds.timestep.setNt(Nt)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=1e20)


# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
 
temperature = T_initial * np.ones((len(times), len(radius)))
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)

# Hot-tail grid settings
# Set initial Maxwellian @ T = 1 keV, n = 5e19, uniform in radius
if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=0, T0=T_initial)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)


# Use the linear solver
#ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Use the new nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setTolerance(reltol=0.01)
ds.solver.setMaxIterations(maxiter = 100)
ds.solver.setVerbose(True)


ds.other.include('fluid', 'lnLambda','nu_s','nu_D')


# Save settings to HDF5 file
ds.save('dream_settings.h5')


###########
# RESTART #
###########

# Duration of restarted self-consistent simulation:
Tmax = 1e-5
Nt = 10

ds2 = DREAMSettings(ds)

#ds2.fromOutput('output.h5', ignore=['n_i'])
ds2.fromOutput('output.h5')

ds2.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds2.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = E_wall*6.28318530718, wall_radius=radius_wall)
ds2.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

ds2.timestep.setTmax(Tmax)
ds2.timestep.setNt(Nt)

ds2.save('dream_settings_restart.h5')

