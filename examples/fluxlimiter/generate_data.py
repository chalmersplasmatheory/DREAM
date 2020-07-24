#!/usr/bin/env python3
#
# This example sets up and runs convergence scans for the basic runaway example
# using three different interpolation schemes for advection coefficients:
#   CENTRED: linear interpolation (second-order accurate, prone to oscillation)
#   QUICK: quadratic upwind interpolation (third-order accurate)
#   MUSCL: a non-linear flux limited method (first- to second-order accurate) 
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
np0 = 1000
nxi0 = 30
ds.hottailgrid.setNxi(nxi0)
ds.hottailgrid.setNp(np0)
ds.hottailgrid.setPmax(pmax)

ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
#ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_PHI_CONST)
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)



# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(1)

# Use the linear solver
ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.solver.setTiming(False)

ds.other.include('fluid/runawayRate')

# Set time stepper
ds.timestep.setTmax(1e-3)
ds.timestep.setNt(10)

#############################
# Set up convergence scan for the three different advection interpolation schemes

npvec = [0,90,130,180,250,350,500,800,np0]
nxivec = [0,3,5,8,14,22,nxi0]
def pIndFunc(index: int, settings: DREAMSettings, baseline):
    #val = baseline + 70*index
    val = npvec[index]
    settings.hottailgrid.setNp(val)
    return settings, val

def xiIndFunc(index: int, settings: DREAMSettings, baseline):
    #val = baseline + 8*index
    val = nxivec[index]
    settings.hottailgrid.setNxi(val)
    return settings, val

def addPScan(cs: ConvergenceScan):
#    cs.addScanParameter(name='hottailgrid.pgrid.np',   baselineValue=np0,  f=pIndFunc,  nvalues=8, startindex=-7)
    cs.addScanParameter(name='hottailgrid.pgrid.np',   baselineValue=np0,  f=pIndFunc,  nvalues=8, startindex=1)
def addXiScan(cs: ConvergenceScan):
#    cs.addScanParameter(name='hottailgrid.xigrid.nxi', baselineValue=nxi0, f=xiIndFunc, nvalues=8, startindex=-7)
    cs.addScanParameter(name='hottailgrid.xigrid.nxi', baselineValue=nxi0, f=xiIndFunc, nvalues=6, startindex=1)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_CENTRED)
ds_CENTRED = DREAMSettings(ds)
cs_CENTRED = ConvergenceScan(ds_CENTRED,outparams=['other.fluid.runawayRate'])
addPScan(cs_CENTRED)
addXiScan(cs_CENTRED)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_QUICK)
ds_QUICK   = DREAMSettings(ds)
cs_QUICK   = ConvergenceScan(ds_QUICK, outparams=['other.fluid.runawayRate'])
addPScan(cs_QUICK)
addXiScan(cs_QUICK)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_SMART)
ds_SMART   = DREAMSettings(ds)
cs_SMART   = ConvergenceScan(ds_SMART, outparams=['other.fluid.runawayRate'])
addPScan(cs_SMART)
addXiScan(cs_SMART)

ds.solver.setType(Solver.NONLINEAR)
ds.solver.setTolerance(1e-2)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_MUSCL)
ds_MUSCL   = DREAMSettings(ds)
cs_MUSCL   = ConvergenceScan(ds_MUSCL, outparams=['other.fluid.runawayRate'])
addPScan(cs_MUSCL)
addXiScan(cs_MUSCL)

'''
ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_SMART_PE)
ds_SMART_PE   = DREAMSettings(ds)
cs_SMART_PE   = ConvergenceScan(ds_SMART_PE, outparams=['other.fluid.runawayRate'])
addPScan(cs_SMART_PE)
addXiScan(cs_SMART_PE)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_MUSCL_PE)
ds_MUSCL_PE   = DREAMSettings(ds)
cs_MUSCL_PE   = ConvergenceScan(ds_MUSCL_PE, outparams=['other.fluid.runawayRate'])
addPScan(cs_MUSCL_PE)
addXiScan(cs_MUSCL_PE)
'''



##############################
# 2. Run convergence scan
##############################

cs_CENTRED.run()
cs_CENTRED.save('convergence_CENTRED.h5')

cs_QUICK.run()
cs_QUICK.save('convergence_QUICK.h5')

cs_SMART.run()
cs_SMART.save('convergence_SMART.h5')

cs_MUSCL.run()
cs_MUSCL.save('convergence_MUSCL.h5')
'''
cs_SMART_PE.run()
cs_SMART_PE.save('convergence_SMART_PE.h5')

cs_MUSCL_PE.run()
cs_MUSCL_PE.save('convergence_MUSCL_PE.h5')
'''


