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

from DREAM import DREAMSettings, runiface
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.OhmicCurrent as JOhm
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.TransportSettings as Transport

ds = DREAMSettings()
#ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

# Physical parameters
E = 0.09       # Electric field strength (V/m)
n = 5e18    # Electron density (m^-3)
T = 1000     # Temperature (eV)
Ip = 700e3  # Plasma current (A)

# Grid parameters
pMax = 1    # maximum momentum in units of m_e*c
Np   = 300  # number of momentum grid points
Nxi  = 20   # number of pitch grid points
tMax = 5e-1 # simulation time in seconds
Nt   = 200   # number of time steps

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

# Hot-tail grid settings
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)

# Set boundary condition type at pMax
#ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)
#ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_ZERO)

# Enable frozen current mode
ds.eqsys.f_hot.transport.setFrozenCurrentMode(Transport.FROZEN_CURRENT_MODE_BETAPAR, Ip_presc=Ip)
#ds.eqsys.f_hot.transport.prescribeDiffusion(1)
ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)
ds.eqsys.n_re.transport.setFrozenCurrentMode(Transport.FROZEN_CURRENT_MODE_BETAPAR, Ip_presc=Ip)
#ds.eqsys.n_re.transport.prescribeDiffusion(1)
ds.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
#ds.eqsys.j_ohm.setCorrectedConductivity(JOhm.CORRECTED_CONDUCTIVITY_DISABLED)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setWallRadius(0.22)
ds.radialgrid.setNr(10)

# Set solver type
ds.solver.setType(Solver.NONLINEAR) # semi-implicit time stepping
#ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.preconditioner.setEnabled(False)
ds.solver.setVerbose(True)
ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)

# include otherquantities to save to output
ds.other.include('fluid','scalar')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

runiface(ds, 'output.h5', quiet=False)

