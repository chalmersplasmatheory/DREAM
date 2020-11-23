#!/usr/bin/env python3
#
# This example first sets up a simple runaway scenario, which is
# then passed to a DREAM.ConvergenceScan object. The convergence
# scan is configured to apply to the most relevant resolution
# parameters for this scenario. We will use the runaway rate as
# a measure of convergence, and we will consider it to be converged
# when it no longer varies significantly. Once set up, the convergence
# scan is also executed and its result are presented.
#
# Run as
#
#   $ ./run.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.ConvergenceScan import ConvergenceScan
from DREAM.ConvergenceScanPlot import ConvergenceScanPlot
from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions


###############################
# 1. Set up baseline scenario
###############################
ds = DREAMSettings()

#E = 0.3     # Electric field strength (V/m)
E = 6.745459970079014
n = 5e19    # Electron density (m^-3)
#T = 1e3     # Temperature (eV)
T = 100

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.avalanche = False

# Hot-tail grid settings
pmax = 2
ds.hottailgrid.setNxi(15)
ds.hottailgrid.setNp(300)
ds.hottailgrid.setPmax(pmax)

ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST)
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_QUICK)


# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setWallRadius(0.22)
ds.radialgrid.setNr(1)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.output.setTiming(stdout=True)

ds.other.include('fluid/runawayRate')

# Set time stepper
ds.timestep.setTmax(1e-3)
ds.timestep.setNt(10)

##############################
# 2. Set up convergence scan
##############################
cs = ConvergenceScan(ds, inparams=['hottailgrid.np', 'hottailgrid.nxi', 'nt'], outparams=['other.fluid.runawayRate'])

cs.run()
cs.save('convergence.h5')

#####################################
# 3. Do stuff with the finished scan
#####################################
csp = ConvergenceScanPlot(cs)
csp.plot(normalized=True)

# Alternatively, we could read back the saved file
# and plot it
#csp = ConvergenceScanPlot('convergence.h5')
#csp.plot(normalized=True)

