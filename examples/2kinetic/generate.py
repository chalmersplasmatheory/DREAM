#!/usr/bin/env python3
#
# This example illustrates how to run DREAM with both the hot and
# runaway electron grids enabled simultaneously (i.e. two kinetic
# grids).
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
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaway
import DREAM.Settings.Solver as Solver


ds = DREAMSettings()

E = 1       # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 1e3     # Temperature (eV)

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.avalanche = Runaway.AVALANCHE_MODE_NEGLECT

nxi = 50
np  = 400
# Hot-tail grid settings
pmaxHot = 2
ds.hottailgrid.setNxi(nxi)
ds.hottailgrid.setNp(np)
ds.hottailgrid.setPmax(pmaxHot)

# Runaway grid settings
pmaxRE = 2*pmaxHot
ds.runawaygrid.setNxi(1*nxi)
ds.runawaygrid.setNp(2*np)
ds.runawaygrid.setPmax(pmaxRE)
#ds.runawaygrid.setEnabled(False)

ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST)
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_DPHI_CONST)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(1)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)
#ds.solver.setType(Solver.NONLINEAR)

ds.output.setFilename('output.h5')
ds.output.setTiming(True)

ds.other.include('fluid/runawayRate')

# Set time stepper
ds.timestep.setTmax(1e-3)
ds.timestep.setNt(10)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

