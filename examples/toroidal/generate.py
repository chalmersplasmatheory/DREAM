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
import scipy.constants
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.RadialGrid as RGrid
import DREAM.Settings.Solver as Solver

ds = DREAMSettings()

# Physical parameters
E = 3       # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)

# Grid parameters
pMax = 0.4    # maximum momentum in units of m_e*c
Np   = 300  # number of momentum grid points
Nxi  = 60   # number of pitch grid points
nr   = 3    # number of radial grid points
tMax = 2e-3 # simulation time in seconds
Nt   = 10   # number of time steps

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
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0) # F=0 outside the boundary

ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_EXPLICIT)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid

ds.radialgrid.setType(RGrid.TYPE_ANALYTIC_TOROIDAL)

mu0 = scipy.constants.mu_0
R0 = 0.68
#R0 = np.inf
a = 0.22
rref = np.linspace(0, a, 20)
IpRef = 1e6 # reference plasma current which generates the poloidal magnetic field (assumed uniform profile)
psiref = -mu0 * IpRef * (1-(rref/a)**2) * a

rDelta = np.linspace(0, a, 20)
Delta  = np.linspace(0, 0.1*a, rDelta.size)
rdelta = np.linspace(0, a, 20)
delta  = np.linspace(0, 0.2, rdelta.size)
ds.radialgrid.setShaping(psi=psiref, rpsi=rref, GOverR0=5.0, kappa=1.5, delta=delta, rdelta=rdelta, Delta=Delta, rDelta=rDelta)
#ds.radialgrid.setShaping(psi=0,G=5.0, kappa=1.5, delta=0.2, Delta=Delta, rDelta=rDelta)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(a)
ds.radialgrid.setMajorRadius(R0)
ds.radialgrid.setNr(nr)

ds.radialgrid.visualize(nr=8, ntheta=40)

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

