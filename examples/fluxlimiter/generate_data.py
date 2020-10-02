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
import DREAM.Settings.Equations.RunawayElectrons as Runaways
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
ds.eqsys.n_re.avalanche = Runaways.AVALANCHE_MODE_NEGLECT

# Hot-tail grid settings
pmax = .3
np0 = 500
nxi0 = 40
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
ds.output.setTiming(stdout=False, file=True)

ds.other.include('fluid/runawayRate')

# Set time stepper
ds.timestep.setTmax(1e-3)
ds.timestep.setNt(10)

#############################
# Set up convergence scan for the different advection interpolation schemes

npvec = [0,20,30,40,50,60,70,80,90,100,110,120,130,150,200,350,np0]
nxivec = [0,3,5,8,15,25,nxi0]
def pIndFunc(index: int, settings: DREAMSettings, baseline):
    val = npvec[index]
    settings.hottailgrid.setNp(val)
    return settings, val

def xiIndFunc(index: int, settings: DREAMSettings, baseline):
    val = nxivec[index]
    settings.hottailgrid.setNxi(val)
    return settings, val

def addPScan(cs: ConvergenceScan):
    cs.addScanParameter(name='hottailgrid.pgrid.np',   baselineValue=np0,  f=pIndFunc,  nvalues=16, startindex=1)
def addXiScan(cs: ConvergenceScan):
    cs.addScanParameter(name='hottailgrid.xigrid.nxi', baselineValue=nxi0, f=xiIndFunc, nvalues=6, startindex=1)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_CENTRED)
ds_CENTRED = DREAMSettings(ds, chain=False)
cs_CENTRED = ConvergenceScan(ds_CENTRED,outparams=['other.fluid.runawayRate'])
addPScan(cs_CENTRED)
addXiScan(cs_CENTRED)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_QUICK)
ds_QUICK   = DREAMSettings(ds, chain=False)
cs_QUICK   = ConvergenceScan(ds_QUICK, outparams=['other.fluid.runawayRate'])
addPScan(cs_QUICK)
addXiScan(cs_QUICK)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_UPWIND)
ds_UPWIND   = DREAMSettings(ds, chain=False)
cs_UPWIND   = ConvergenceScan(ds_UPWIND, outparams=['other.fluid.runawayRate'])
addPScan(cs_UPWIND)
addXiScan(cs_UPWIND)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_UPWIND_2ND_ORDER)
ds_UPWIND2   = DREAMSettings(ds, chain=False)
cs_UPWIND2   = ConvergenceScan(ds_UPWIND2, outparams=['other.fluid.runawayRate'])
addPScan(cs_UPWIND2)
addXiScan(cs_UPWIND2)


ds.solver.setType(Solver.NONLINEAR)
ds.solver.setTolerance(1e-2)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_SMART)
ds_SMART   = DREAMSettings(ds, chain=False)
cs_SMART   = ConvergenceScan(ds_SMART, outparams=['other.fluid.runawayRate'])
addPScan(cs_SMART)
addXiScan(cs_SMART)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_MUSCL)
ds_MUSCL   = DREAMSettings(ds, chain=False)
cs_MUSCL   = ConvergenceScan(ds_MUSCL, outparams=['other.fluid.runawayRate'])
addPScan(cs_MUSCL)
addXiScan(cs_MUSCL)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_OSPRE)
ds_OSPRE   = DREAMSettings(ds, chain=False)
cs_OSPRE   = ConvergenceScan(ds_OSPRE, outparams=['other.fluid.runawayRate'])
addPScan(cs_OSPRE)
addXiScan(cs_OSPRE)

ds.eqsys.f_hot.setAdvectionInterpolationMethod(fluxlimiterdamping=1, ad_int=FHot.AD_INTERP_TCDF)
ds_TCDF   = DREAMSettings(ds, chain=False)
cs_TCDF   = ConvergenceScan(ds_TCDF, outparams=['other.fluid.runawayRate'])
addPScan(cs_TCDF)
addXiScan(cs_TCDF)

##############################
# 2. Run convergence scan
##############################
'''
cs_CENTRED.run()
cs_CENTRED.save('convergence_CENTRED.h5')

cs_QUICK.run()
cs_QUICK.save('convergence_QUICK.h5')

cs_UPWIND.run()
cs_UPWIND.save('convergence_UPWIND.h5')

cs_UPWIND2.run()
cs_UPWIND2.save('convergence_UPWIND2.h5')

cs_SMART.run()
cs_SMART.save('convergence_SMART.h5')
'''
cs_MUSCL.run()
cs_MUSCL.save('convergence_MUSCL.h5')

cs_OSPRE.run()
cs_OSPRE.save('convergence_OSPRE.h5')

cs_TCDF.run()
cs_TCDF.save('convergence_TCDF.h5')

