#!/usr/bin/env python3
#
# This example shows how to use the ionization-based adaptive time
# stepper. The example sets up a simple 0D plasma with argon
# impurities injected at t=0. The impurities will ionize and lead to
# a rapid cooling of the plasma. Initially, the time scales should
# be very fast but gradually become longer as the simulation
# proceeds, serving as a good example of the purpose of the
# ionization-based adaptive time stepper.
#
# Run as
#
#   $ ./generate.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM import DREAMSettings, runiface
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Solver as Solver


ds = DREAMSettings()

#############################
# Set simulation parameters #
#############################
T_initial = 500     # initial temperature (eV)
Ip_initial = 400e5  # initial plasma current (A)

TMAX = 1e-3         # simulation time (s)
NR = 1              # number of radial grid points

B0 = 5              # magnetic field strength (T)
a = 1.0             # plasma minor radius (m)
b = 1.1             # wall minor radius (m)
R0 = 3.0            # plasma major radius (m)

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(b)
ds.radialgrid.setNr(NR)

# Set up time stepper
ds.timestep.setIonization(dt0=1e-7, dtmax=4e-6, tmax=TMAX)
ds.timestep.setMinSaveTimestep(3e-7)

# Add ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=1e19)
ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=2e18)

# Set initial plasma current/electric field
ds.eqsys.j_ohm.setInitialProfile(1, Ip0=Ip_initial)
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0, R0=R0)
# Set initial temperature
ds.eqsys.T_cold.setType(T_cold.TYPE_SELFCONSISTENT)
ds.eqsys.T_cold.setInitialProfile(T_initial)

# Set runaway rates
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Disable kinetic grids
ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Use the nonlinear solver
#ds.solver.setType(Solver.NONLINEAR)
#ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_GMRES)
#ds.solver.setTolerance(reltol=0.01)
#ds.solver.setMaxIterations(maxiter = 100)
ds.solver.setVerbose(False)

ds.other.include('fluid')

# Run the simulation
ds.save('settings.h5')
runiface(ds, 'output.h5')


