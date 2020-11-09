#!/usr/bin/env python3
#
# This example shows how to run a combined fluid-kinetic simulation with
# with both the hot-tail and runaway electron grids.
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
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.RunawayElectrons as Runaways


ds = DREAMSettings()

E = 0.6     # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 1e3     # Temperature (eV)

re_enabled = True

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable hot-tail grid
ds.hottailgrid.setEnabled(False)

ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_ULTRA_RELATIVISTIC
#ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL

# Set initial hot electron Maxwellian
#ds.eqsys.f_hot.setInitialProfiles(n0=2*n, T0=T)

# Include Dreicer and avalanche
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Disable runaway grid
pmax_re = 0.5
if re_enabled:
    ds.runawaygrid.setNxi(50)
    ds.runawaygrid.setNp(100)
    ds.runawaygrid.setPmax(pmax_re)

    # Use flux limiters
    ds.eqsys.f_re.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)

else:
    ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(1)

# Use the linear solver
#ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setVerbose(True)

ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(1e-1)
ds.timestep.setNt(20)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

