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

sys.path.append('/home/qhli/DREAM/DREAM/py')
try:
    import DREAM
except ModuleNotFoundError:
    sys.path.append('/home/qhli/DREAM/DREAM/py')
    import DREAM

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaway
import DREAM.Settings.RadialGrid as RGrid
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.DistributionFunction as DistFunc

# 配置输出文件名称
output_base_name = 'output15'  # 修改这一行就能改变所有输出文件名                                                                                                                                                                  
output_settings_file = f'dream_settings_{output_base_name}.h5'
output_result_file = f'{output_base_name}.h5'

ds = DREAMSettings()

E = 16       # Electric field strength (V/m)
n = 2e19    # Electron density (m^-3)
T = 5e2     # Temperature (eV)

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
#ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)
ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED, Z0=12, n=n)
#ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
#ds.eqsys.n_re.avalanche = Runaway.AVALANCHE_MODE_NEGLECT
ds.eqsys.n_re.setDreicer(Runaway.DREICER_RATE_DISABLED)
ds.eqsys.n_re.setAvalanche(Runaway.AVALANCHE_MODE_KINETIC, pCutAvalanche = 0.01)

ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED 
ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS 

nxi = 20
nP  = 100
# Hot-tail grid settings
pmaxHot = 0.94
ds.hottailgrid.setNxi(nxi)
ds.hottailgrid.setNp(nP)
ds.hottailgrid.setPmax(pmaxHot)

# Runaway grid settings
pmaxRE = 50*pmaxHot
ds.runawaygrid.setNxi(1*nxi)
ds.runawaygrid.setNp(2*nP)
ds.runawaygrid.setPmax(pmaxRE)

ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)
ds.runawaygrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)


# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)
ds.eqsys.f_hot.setParticleSource(particleSource=FHot.PARTICLE_SOURCE_IMPLICIT)
#ds.eqsys.f_re.setBoundaryCondition(DistFunc.BC_PHI_CONST)
ds.eqsys.f_re.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)
#ds.eqsys.n_re.setEceff(Runaway.COLLQTY_ECEFF_MODE_EC_TOT)


a=0.9                                                                                                                                                                                                                                                                                                 
b=0.99
R0=1.0
B0=2.5
nr=6

ds.radialgrid.setType(RGrid.TYPE_ANALYTIC_TOROIDAL)
ds.radialgrid.setMajorRadius(R0)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setShaping(psi=0.0001, GOverR0=B0)
#ds.radialgrid.setNr(5)
r_f = np.linspace(0.6, a, nr+1)   
r_f *= a/r_f[-1]
ds.radialgrid.setCustomGridPoints(r_f)
print(r_f)
ds.radialgrid.setWallRadius(b)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)
#ds.solver.setType(Solver.NONLINEAR)
#ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
ds.solver.setVerbose(False)

ds.output.setFilename(output_result_file)
ds.output.setTiming(True)

ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(5e-5)
ds.timestep.setNt(50)
##ds.timestep.setNumberOfSaveSteps(50)

# Save settings to HDF5 file
ds.save(output_settings_file)
do = DREAM.runiface(ds, output_result_file, quiet=False)