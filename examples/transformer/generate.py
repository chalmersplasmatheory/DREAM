#!/usr/bin/env python3
#
# This example illustrates the use of the 'TRANSFORMER' boundary
# condition for the poloidal flux equation, whereby a resistive wall
# is assumed, along with a non-zero applied loop voltage via the
# transformer (passing through major radius R = 0).
#
# Run as
#
#   $ ./generate.py
#   $ ../../build/iface/dreami dream_settings.h5
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions

ds = DREAMSettings()

# Physical parameters
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)

# Grid parameters
tMax = 5e-1 # simulation time in seconds
Nt   = 40   # number of time steps
Nr   = 30

# Tokamak parameters
a    = 0.9          # (m)
b    = 1.2          # (m)
B0   = 3.1          # (T)
R0   = 2.96         # (m)
tau_wall = 0.01     # (s)

# Set E_field
#ds.eqsys.E_field.setInitialProfile(0)
ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=0.4/R0, inverse_wall_time=1/tau_wall, R0=R0)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable kinetic grids
ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(b)
ds.radialgrid.setNr(Nr)

# Set solver type
ds.solver.setType(Solver.NONLINEAR)

# Include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

