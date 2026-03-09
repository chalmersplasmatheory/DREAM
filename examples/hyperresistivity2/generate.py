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
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions

from DREAM import DREAMSettings, DREAMOutput, DREAMException, runiface

ds = DREAMSettings()

# Physical parameters
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)
a = 0.22

tMax = 5e-6 # simulation time in seconds
Nt   = 100   # number of time steps

# Set E_field
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(
    bctype=EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0,
    R0=0.88
)

# Set current density
r = np.linspace(0, a)
j0 = 1e6 * (1-(r/a)**2)
ds.eqsys.j_ohm.setInitialProfile(j0, r, Ip0=250e3)

ds.eqsys.psi_p.setHyperresistivity([1e-4], radius=[0], times=[0])

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setDreicer(True)
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_FLUID)
ds.eqsys.n_re.setInitialProfile(1e5)

# Disable kinetic grids
ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(a*1.1)
ds.radialgrid.setNr(20)

# Set solver type
ds.solver.setType(Solver.NONLINEAR) # semi-implicit time stepping

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

do = runiface(ds, 'output.h5')
