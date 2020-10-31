#!/usr/bin/env python3
#
# This script generates a simple CODE-like Dreicer runaway simulation
# with additional pitch scattering due to a magnetic ripple included.
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
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

ds = DREAMSettings()

# Physical parameters
E = .1       # Electric field strength (V/m)
n = 1e19    # Electron density (m^-3)
T = 1000     # Temperature (eV)

# Grid parameters
pMax = 12    # maximum momentum in units of m_e*c
#Np   = 400  # number of momentum grid points
#Nxi  = 30   # number of pitch grid points
Np   = 200  # number of momentum grid points
Nxi  = 20   # number of pitch grid points
tMax = 2e-1 # simulation time in seconds
Nt   = 20   # number of time steps

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
ds.hottailgrid.setBiuniformGrid(psep=0.5, npsep=50)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)

# Set boundary condition type at pMax
#ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary

# Magnetic ripple
ds.eqsys.f_hot.setRipple(deltacoils=0.35, m=1, n=10, dB_B=1e-4)
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(3)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(1)

# Set solver type
ds.solver.setType(Solver.NONLINEAR) # semi-implicit time stepping
ds.solver.tolerance.set(reltol=1e-4)
ds.solver.setVerbose(False)

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

