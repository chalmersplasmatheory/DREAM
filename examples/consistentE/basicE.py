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
import DREAM.Settings.Equations.HotElectronDistribution as FHot
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

# time resolution of restarted simulation
Tmax_restart = 1e-3 # simulation time in seconds
Nt_restart = 20     # number of time steps

B0 = 5              # magnetic field strength in Tesla
E_initial = 60      # initial electric field in V/m
E_wall = 0.0        # boundary electric field in V/m
T_initial = 6       # initial temperature in eV

Tmax_init = 1e-5    # simulation time in seconds
Nt_init = 3         # number of time steps
Nr = 4              # number of radial grid points
Np = 200            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 0.03         # maximum momentum in m_e*c
times  = [0]        # times at which parameters are given
radius = [0, 1]     # span of the radial grid
radius_wall = 1.5   # location of the wall 

T_selfconsistent    = True
hotTailGrid_enabled = True

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)
# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=1e20)
ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e20)
#ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
#ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=1e20)


# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
ds.eqsys.E_field.setBoundaryCondition(wall_radius=radius_wall)

#if T_selfconsistent:
#    temperature = T_initial * np.ones(len(radius))
#    ds.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)
#    ds.eqsys.T_cold.setInitialProfile(temperature=temperature,radius=radius)
#else:
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
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)


# Disable runaway grid
ds.runawaygrid.setEnabled(False)


# Use the linear solver
#ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Use the new nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
#ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_GMRES)
ds.solver.setTolerance(reltol=0.01)
ds.solver.setMaxIterations(maxiter = 100)
ds.solver.setVerbose(True)


ds.other.include('fluid', 'lnLambda','nu_s','nu_D')

ds.output.setFilename('output_init.h5')

# Save settings to HDF5 file
ds.save('init_settings.h5')


###########
# RESTART #
###########

ds2 = DREAMSettings(ds)

ds2.fromOutput('output_init.h5')

ds2.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds2.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = E_wall*2*np.pi, wall_radius=radius_wall)
if T_selfconsistent:
    ds2.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

ds2.timestep.setTmax(Tmax_restart)
ds2.timestep.setNt(Nt_restart)

# I need to do the below for the boundary condition to be changed
#ds2.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)

ds2.save('restart_settings.h5')

