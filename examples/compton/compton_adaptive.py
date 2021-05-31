#!/usr/bin/env python3
#
# This example shows how to set up a self-consistent fluid DREAM run,
# where no kinetic equations are solved, but the electric field and
# temperature are evolved self-consistently.
#
# Run as
#
#   $ ./compton.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
from DREAM.DREAMOutput import DREAMOutput
from DREAM import runiface
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.RunawayElectrons as RE
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.TimeStepper as TimeStep


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

run_ioniz=True
run_expdecay=True
run_CQ=True

filename_ending="Tf100-50_Nr101"

ds = DREAMSettings()

# set collision settings
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
#ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_NEGLECT
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
#ds.collisions.lnlambda = Collisions.LNLAMBDA_CONSTANT
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL

#############################
# Set simulation parameters #
#############################
Tmax_PostTQ = 40e-3

Tmax_restart_expdecay = 6e-3
Nt_restart_expdecay = 30

# time resolution of restarted simulation
Tmax_restart = 0.3e-6 # simulation time in seconds
Nt_restart = 100     # number of time steps

n_D = 1e20
n_D_inj = 40*n_D
n_Z = 0.08*n_D

B0 = 5.3            # magnetic field strength in Tesla
E_initial = 0.00032 # initial electric field in V/m
E_wall = 0.0001        # boundary electric field in V/m
# NOTE: it does not work to have self-consistent E-field with prescribed BC with E_wall=0, 
# since that leads to Psi_wall=0 constantly, which does not work when you have a relative tolerance
T_initial = 20e3    # initial temperature in eV
T_final = 50
t0=1e-3

Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps
Nr = 101             # number of radial grid points
Np = 200            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 1.0          # maximum momentum in m_e*c
times  = [0]        # times at which parameters are given
times_T=np.linspace(0,8.1e-3)
radius = [0, 2]     # span of the radial grid
radialgrid = np.linspace(radius[0],radius[-1],Nr)
radius_wall = 2.15  # location of the wall 

T_selfconsistent    = True
hotTailGrid_enabled = False

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setWallRadius(radius_wall)
ds.radialgrid.setNr(Nr)

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)

# Set ions
#density_D = n_D*np.ones(len(radius))
#density_Ne = n_Ne*np.ones(len(radius))
density_D = n_D
density_D_inj = n_D_inj
density_Z = n_Z

ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=density_D)
ds.eqsys.n_i.addIon(name='D_inj', Z=1, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=density_D_inj)
ds.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=density_Z)
#ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
#ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=1e20)


# Set E_field 
"""
do=DREAMOutput('output_init.h5')
conductivity=do.other.fluid.conductivity.getData()
jprof=(1-(1-0.001**(1/0.41))*(radialgrid/radialgrid[-1])**2)**0.41
efield=1.81e6*jprof/conductivity[-1,:]
"""

efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)

# Set runaway generation rates
ds.eqsys.n_re.setCompton(RE.COMPTON_RATE_ITER_DMS)
ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)
ds.eqsys.n_re.setEceff(RE.COLLQTY_ECEFF_MODE_CYLINDRICAL)

# temperature = T_initial * np.ones((len(times), len(radius)))
#temperature = T_final+(T_initial - T_final) * np.exp(-times_T/t0).reshape(-1,1) * np.ones((len(times_T), len(radius)))
temp_prof=(1-0.99*(radialgrid/radialgrid[-1])**2).reshape(1,-1)
temperature = T_final+(T_initial*temp_prof - T_final) * np.exp(-times_T/t0).reshape(-1,1)
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times_T, radius=radialgrid)

if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

#nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
#ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=0, T0=T_initial)
#ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)


# Use the linear solver
# ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Use the nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_LU)
ds.solver.tolerance.set(reltol=1e-4)
ds.solver.setMaxIterations(maxiter = 500)
# ds.solver.setVerbose(True)


ds.other.include('fluid', 'lnLambda','nu_s','nu_D')


# Save settings to HDF5 file
ds.save('init_settings_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5')
runiface(ds, 'output_init_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5', quiet=False)


#######################
# RESTART set current #
#######################

do=DREAMOutput('output_init_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5')
conductivity=do.other.fluid.conductivity.getData()
jprof=(1-(1-0.001**(1/0.41))*(radialgrid/radialgrid[-1])**2)**0.41
# efield=1.81e6*jprof/conductivity[-1,:]
efield=1.69e6*jprof/conductivity[-1,:]

ds.eqsys.E_field.setPrescribedData(efield=efield, radius=radialgrid)

# Save settings to HDF5 file
ds.save('init_settings_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5')
runiface(ds, 'output_init_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5', quiet=False)


#################
# RESTART ioniz #
#################

ds2 = DREAMSettings(ds)


ds2.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds2.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall_R0 = E_wall*2*np.pi)

ds2.timestep.setTmax(Tmax_restart)
ds2.timestep.setNt(Nt_restart)

ds2.save('ioniz_restart_settings_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5')
if run_ioniz:
    runiface(ds2, 'output_restart_ioniz_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5', quiet=False)

####################
# RESTART expdecay #
####################
ds3 = DREAMSettings(ds2)

ds3.timestep.setTmax(Tmax_restart_expdecay)
ds3.timestep.setNt(Nt_restart_expdecay)

ds3.save('expdecay_restart_settings_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5')
if run_expdecay:
    runiface(ds3, 'output_restart_expdecay_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5', quiet=True)

###################
# RESTART CQ init #
###################

tMaxTest = 1e-6
dtTest = 1e-7

ds4_bm = DREAMSettings(ds3)
ds4_bm.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)
ds4_bm.timestep.setTmax(tMaxTest)
ds4_bm.timestep.setNt(round(tMaxTest/dtTest))
runiface(ds4_bm, 'output_test1.h5', quiet=False)
ds4_bm.timestep.setNt(10*round(tMaxTest/dtTest))
runiface(ds4_bm, 'output_test10.h5', quiet=False)

ds4 = DREAMSettings(ds3)
if T_selfconsistent:
    ds4.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)
ds4.timestep.setType(TimeStep.TYPE_ADAPTIVE)
ds4.timestep.setTmax(tMaxTest)
#ds4.timestep.setTmax(Tmax_PostTQ)
ds4.timestep.setCheckInterval(1)
ds4.timestep.setRelativeTolerance(1e-1)
ds4.timestep.setVerbose(False)
ds4.timestep.setDt(1e-8)

ds4.save('PostTQ_restart_settings_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5')
if run_CQ:
    runiface(ds4, 'output_restart_PostTQ_init_nNe'+str(n_Z)+'nD_inj'+str(n_D_inj)+filename_ending+'.h5', quiet=False)
