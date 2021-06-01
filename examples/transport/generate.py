#!/usr/bin/env python3
#
# This script is used to test radial transport of n_re with a (constant)  
# scalar diffusion coefficient. Hot-tail RE generated with exp. temperature drop.
#
# By Ida Svenningsson, 2020
#
# ###################################################################

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
from DREAM.DREAMOutput import DREAMOutput
from DREAM import runiface

import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.RunawayElectrons as RE
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

from DREAM import DREAMIO

ds = DREAMSettings()

# set collision settings

######################
# COLLISION SETTINGS #
######################

ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL 
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED # NON_SCREENED, PARTIALLY_SCREENED, COMPLETELY_SCREENED
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL

run_init = True
run_exp = True

#transport_mode = Transport.TRANSPORT_PRESCRIBED
transport_mode = Transport.TRANSPORT_RECHESTER_ROSENBLUTH

#############################
# Set simulation parameters #
#############################

# initial run (to get correct E-field profile)
Tmax_init = 1e-11 
Nt_init = 2       

# Exponential temperature drop
Tfinal_exp = 50 
t0_exp = .5e-3 
Tmax_exp = 10e-3 
Nt_exp = 3000
times_exp = np.linspace(0,Tmax_exp,Nt_exp) 

n_D = 1e20 # originally present in the plasma

B0 = 5              # magnetic field strength in Tesla
E_initial = 5e-3    # initial electric field in V/m
E_wall = 0.0        # boundary electric field in V/m
T_initial = 20e3    # initial temperature in eV
T_in_back = 10      # initial T_cold value (?)
jTot = 1.69e6

Nr_kin  = 15        # number of radial grid points
Np      = 100        # number of momentum grid points
Nxi     = 1         # number of pitch grid points
pMax    = 3         # maximum momentum in m_e*c
times   = [0]       # times at which parameters are given
radius  = [0, 2]    # span of the radial grid
radius_wall = 2.15  # location of the wall 

diffusion_coeff = 100 # m/s^2   -- Diffusion coefficient

hotTailGrid_enabled = True
if hotTailGrid_enabled == False and transport_mode == Transport.TRANSPORT_RECHESTER_ROSENBLUTH:
    print('WARNING: Using Rechester-Rosenbluth transport requires f_hot. Enabling hot-tail grid...')

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setWallRadius(radius_wall)
ds.radialgrid.setNr(Nr_kin)

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)

# Set ions
density_D = n_D

#ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=density_D)
#ds.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=density_Z)
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=density_D)
#ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_PRESCRIBED_NEUTRAL, n=n_Z)


# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)

radialgrid = np.linspace(radius[0],radius[-1],Nr_kin)
temp_prof=(1-0.99*(radialgrid/radialgrid[-1])**2).reshape(1,-1)
temperature_init = Tfinal_exp+(T_initial*temp_prof - Tfinal_exp)

ds.eqsys.T_cold.setPrescribedData(temperature=temperature_init, times=[0], radius=radialgrid)

ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID)
ds.hottailgrid.setEnabled(False) # To be enabled later

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Use the new nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.tolerance.set(reltol=1e-4)
ds.solver.setMaxIterations(maxiter = 100)
ds.solver.setVerbose(False)
ds.output.setTiming(False)
ds.other.include('fluid', 'transport')

if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setEnabled(True)
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)
    nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    #ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial*.99, rT0=0, T0=T_initial)
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial*.99, rT0=radialgrid, T0=temperature_init[0,:])
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)

#########################################
# XXX: Prescribe diffusive transport
#########################################
if transport_mode == Transport.TRANSPORT_PRESCRIBED:
    ds.eqsys.n_re.transport.prescribeDiffusion(drr=diffusion_coeff)
    ds.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
    # with hyperresistivity:
#    ds.eqsys.psi_p.transport.prescribeDiffusion(drr=1e-5)
#    ds.eqsys.psi_p.transport.setBoundaryCondition(Transport.BC_CONSERVATIVE)
elif transport_mode  == Transport.TRANSPORT_RECHESTER_ROSENBLUTH and hotTailGrid_enabled:
    ds.eqsys.f_hot.transport.setMagneticPerturbation(1e-5)
    ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)
#########################################

########
# init # 
########
# To get the right initial current profile

if run_init:
    ds.save('initsim.h5')
    runiface(ds,f'out_1init.h5')

######################################
# RE-SCALE E-FIELD FOR RIGHT CURRENT #
######################################
do=DREAMOutput(f'out_1init.h5')
conductivity=do.other.fluid.conductivity.getData()
jprof=(1-(1-0.001**(1/0.41))*(radialgrid/radialgrid[-1])**2)**0.41
efield=jTot*jprof/conductivity[-1,:]
ds.eqsys.E_field.setPrescribedData(efield=efield, radius=radialgrid)

ds.save(f'settings_1init.h5')        
if run_init:
    runiface(ds,f'out_1init.h5')

#######
# exp # 
#######
# Start expdecay

ds3 = DREAMSettings(ds)
if run_exp:
    ds3.fromOutput(f'out_1init.h5')

temperature_exp = Tfinal_exp+(T_initial*temp_prof - Tfinal_exp) * np.exp(-times_exp/t0_exp).reshape(-1,1)
ds3.eqsys.T_cold.setPrescribedData(temperature=temperature_exp, times=times_exp, radius=radialgrid)
ds3.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds3.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall_R0 = E_wall*2*np.pi)

ds3.timestep.setTmax(Tmax_exp)
ds3.timestep.setNt(Nt_exp)

ds3.save(f'settings_2exp.h5')
if run_exp:
    runiface(ds3,f'out_2exp.h5')
        
