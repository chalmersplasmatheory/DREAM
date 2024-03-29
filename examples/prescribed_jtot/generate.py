#!/usr/bin/env python3
#
# This example shows how to set up a simple case where the current
# density is prescribed. The simulation uses a constant temperature and
# density.
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
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions

ds = DREAMSettings()
#ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

# Physical parameters
E = 6       # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)

# Grid parameters
pMax = 1    # maximum momentum in units of m_e*c
Np   = 300  # number of momentum grid points
Nxi  = 20   # number of pitch grid points
tMax = 1e-3 # simulation time in seconds
Nt   = 2   # number of time steps

with_fhot = False

# Set E_field
#ds.eqsys.E_field.setPrescribedData(E)
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0, R0=0.68)
ds.eqsys.j_ohm.setInitialProfile(1e5, Ip0=2e5)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

# Hot-tail grid settings
if with_fhot:
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
else:
    ds.hottailgrid.setEnabled(False)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setWallRadius(0.22)
ds.radialgrid.setNr(10)

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.preconditioner.setEnabled(False)

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

