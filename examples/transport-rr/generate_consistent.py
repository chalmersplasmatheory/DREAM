#!/usr/bin/env python3
#
# This example shows how to set up a simple CODE-like runaway
# scenario in DREAM. The simulation uses a constant temperature,
# density and electric field, and generates a runaway current
# through the electric field acceleration, demonstrating Dreicer generation.
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

from DREAM import runiface
from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport

ds = DREAMSettings()

# Physical parameters
n = 5e19    # Electron density (m^-3)
T = 10e3    # Temperature (eV)

Ip0 = 1e6  # Initial plasma current (A)

# Grid parameters
pMax = 1    # maximum momentum in units of m_e*c
Np   = 150  # number of momentum grid points
Nxi  = 6    # number of pitch grid points
tMax = 2e-3 # simulation time in seconds
Nt   = 30   # number of time steps
Nr   = 4    # number of radial grid points

minor_radius = 0.22     # m
dBOverB = 3e-3  # Magnetic perturbation strength

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(minor_radius)
ds.radialgrid.setNr(Nr)

# Hot-tail grid settings
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)
# Set boundary condition type at pMax
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0) # F=0 outside the boundary
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST) # extrapolate flux to boundary

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.avalanche = Runaways.AVALANCHE_MODE_NEGLECT

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set time stepper
ds.timestep.setTmax(5e-1)
ds.timestep.setNt(4)

# Set initialization E_field
ds.eqsys.E_field.setPrescribedData(Ip0 / 1.56e8)

# include otherquantities to save to output
ds.other.include('fluid')

runiface(ds, 'init_output.h5', quiet=True)

################################################
######## BEGIN RESTART SIMULATION SETUP ########
################################################

ds_re = DREAMSettings(ds)

# Set time stepper
ds_re.timestep.setTmax(tMax)
ds_re.timestep.setNt(Nt)

ds_re.eqsys.T_cold.setType(T_cold.TYPE_SELFCONSISTENT)
ds_re.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds_re.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = 0, wall_radius=1.1*minor_radius)

# Set Rechester-Rosenbluth transport
# in T_cold
ds_re.eqsys.T_cold.transport.setMagneticPerturbation(dBOverB)
ds_re.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)
# and in f_hot
ds_re.eqsys.f_hot.transport.setMagneticPerturbation(dBOverB)
ds_re.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)

# Set solver type
#ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds_re.solver.setType(Solver.NONLINEAR)
ds_re.solver.setVerbose(True)
ds_re.solver.setLinearSolver(Solver.LINEAR_SOLVER_MUMPS)
ds_re.solver.tolerance.set(reltol=1e-4)
ds_re.solver.setDebug(savejacobian=True, savenumericaljacobian=True, timestep=1,iteration=5)

ds_re.output.setTiming(stdout=True, file=True)
ds_re.output.setFilename('output.h5')

# Save settings to HDF5 file
ds_re.save('dream_settings.h5')
