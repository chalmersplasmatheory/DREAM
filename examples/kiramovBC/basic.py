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
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.TransportSettings as Transport

import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Equations.ElectricField as Efield

from DREAM import DREAMSettings, DREAMOutput, DREAMException, runiface

ds = DREAMSettings()
#ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

# Disable kinetic grids: 
ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Physical parameters
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)

# Grid parameters
tMax = 20e-3 # simulation time in seconds
Nt   = 20# number of time steps

# Set E_field
ds.eqsys.E_field.setPrescribedData(0)


# Set temperature
#x = np.linspace(0.20, 0.50, 30)  # Generate 500 points in the range
#T =100 - x**2 * 10

#ds.eqsys.T_cold.setPrescribedData(temperature=T)
ds.eqsys.T_cold.setType(T_cold.TYPE_SELFCONSISTENT)
ds.eqsys.T_cold.setInitialProfile(T)

ds.eqsys.E_field.setPrescribedData(0.0001)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)


# Set Kiramov Boundary condition for the parallel heat losses 
ds.eqsys.T_cold.transport.prescribeDiffusion(1)
ds.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_KIRAMOV)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.20)
ds.radialgrid.setWallRadius(0.50)
ds.radialgrid.setNr(30)

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.preconditioner.setEnabled(False)

# include otherquantities to save to output
ds.other.include('fluid','scalar')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('/home/votta/DREAM/examples/kiramovBC/output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

do = runiface(ds, 'test.h5')
