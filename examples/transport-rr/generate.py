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

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport

ds = DREAMSettings()

# Physical parameters
E = 6       # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)

# Grid parameters
pMax = 1    # maximum momentum in units of m_e*c
Np   = 300  # number of momentum grid points
Nxi  = 10   # number of pitch grid points
tMax = 2e-2 # simulation time in seconds
Nt   = 20   # number of time steps
Nr   = 10   # number of radial grid points

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.avalanche = Runaways.AVALANCHE_MODE_NEGLECT

# Hot-tail grid settings
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)

# Set boundary condition type at pMax
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0) # F=0 outside the boundary

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(Nr)

# Set Rechester-Rosenbluth transport
ds.eqsys.f_hot.transport.prescribeDiffusion(1e-1)
#ds.eqsys.f_hot.transport.setMagneticPerturbation(1e-4)
ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)
#ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_CONSERVATIVE)

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
#ds.solver.setType(Solver.NONLINEAR)
#ds.solver.setVerbose(True)
#ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MUMPS)

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

