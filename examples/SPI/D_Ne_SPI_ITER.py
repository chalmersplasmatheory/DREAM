#!/usr/bin/env python3
#
# This script shows how to set up an D/Ne SPI scenario in an ITER-like setting.
# The injection can either be separated into one stage with pure D and one stage with pure NE,
# or be made as a single stage injection with a similar total amount of particles.
#
################################################################################################ 

import numpy as np
import scipy as sp
import sys

from scipy import integrate
from scipy.special import kn

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
import DREAM.Settings.Equations.SPI as SPI
import DREAM.Settings.TransportSettings as Transport


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

# Makes sure this script generates exactly the same result every time the same settings are used
np.random.seed(1)

ds = DREAMSettings()


#######################################
# Set numerical simulation parameters #
#######################################

# Choose which parts of the disruption should be simulated
run_init=True # Includes a single-step run to extract the conductivity and 
#                another one to set the efield according to the wanted current profiel

run_injection_init=True # Accelerate the distribution function to carry the right current

run_injection=True # D or D/Ne injection
run_CQ=True # Second Ne injection (if any), beginning of the CQ

# Specify number of restarts to do during the CQ
nCQ_restart_start=2 # Number of CQ restart to start from
nCQ_restart=0 # How many CQ restarts to run

# Temperature and electron distribution settings
T_selfconsistent    = True
hotTailGrid_enabled = False
use_fluid_runaways = False

use_heat_transport=True
use_f_hot_transport=False
transport_CQ_only=False
dBOverB=1e-3/np.sqrt(2)

# Time steps during the various restarts
Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps

Tmax_init2 = 3e-3   
Nt_init2 = 300      

Tmax_injection = 3.4e-3
Nt_injection = 340   

# For single stage
# Tmax_injection = 6e-3
# Nt_injection = 3000  

Tmax_CQ = 17e-3
Nt_CQ = 3000
Tmax_CQ_restart = 0.2e-3
Nt_CQ_restart = 100

# Grid parameters
Nr = 11             # number of radial grid points
Np = 120            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 3            # maximum momentum in m_e*c
times  = [0]        # times at which parameters are given
radius = [0, 2]     # span of the radial grid
dr=(radius[1]-radius[0])/(Nr+1)
radialgrid = np.linspace(radius[0]+dr/2,radius[-1]-dr/2,Nr)
radius_wall = 2.15  # location of the wall 

B0 = 5.3            # magnetic field strength in Tesla

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)
ds.radialgrid.setWallRadius(radius_wall)

#######################################
# Set physical simulation parameters  #
#######################################

# Set E_field 
E_initial = 0.00032 # initial electric field in V/m
E_wall = 0.0        # boundary electric field in V/m
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
ds.eqsys.E_field.setBoundaryCondition()

# Set runaway generation rates
if use_fluid_runaways:
	ds.eqsys.n_re.setCompton(RE.COMPTON_RATE_ITER_DMS)
	ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)

# Set temperature profile
T_initial = 20e3    # initial temperature in eV
temp_prof=(1-0.99*(radialgrid/radialgrid[-1])**2).reshape(1,-1)
temperature = T_initial*temp_prof
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radialgrid)

# Settings for the first SPI (presumably consisting mostly of deuterium)
nShardD=1742 # Number of shards
NinjD=2e24 # Number of atoms
alpha_maxD=0.17 # Divergence angle
abs_vp_meanD=800 # Mean shard speed
abs_vp_diffD=0.2*abs_vp_meanD # Width of the uniform shard speed distribution
molarFractionNe=0 # Molar fraction of neon (the rest is deuterium)

# The shard velocities are set to zero for now, 
# and will be changed later when the injections are supposed to start
if molarFractionNe>0:
	ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardD, Ninj=NinjD, Zs=[1,10], isotopes=[2,0], molarFractions=[1-molarFractionNe,molarFractionNe], ionNames=['D_inj_mix','Ne_inj_mix'], abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxD, shatterPoint=np.array([radius_wall,0,0]))
else:
	ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardD, Ninj=NinjD, Zs=[1], isotopes=[2], molarFractions=[1], ionNames=['D_inj'], abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxD, shatterPoint=np.array([radius_wall,0,0]))

# Settings for the second Neon SPI
nShardNe=50
NinjNe=1e24
alpha_maxNe=0.17
abs_vp_meanNe=200
abs_vp_diffNe=0.2*abs_vp_meanNe

if nShardNe>0:
	ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardNe, Ninj=NinjNe, Zs=[10], isotopes=[0], molarFractions=[1], ionNames=['Ne_inj'], abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxNe, shatterPoint=np.array([2.15,0,0]))

# Set geometrical parameters used to rescale VpVol when calculting the size of the flux surfaces
R=6.2
kappa=1
ds.eqsys.spi.setVpVolNormFactor(R*kappa)

# Set models for advancing the shard motion and ablation
ds.eqsys.spi.setVelocity(SPI.VELOCITY_MODE_PRESCRIBED) # Constant prescribed velocities
ds.eqsys.spi.setAblation(SPI.ABLATION_MODE_FLUID_NGS) # Parks NGS formula based on T_cold
# ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL_GAUSSIAN) # Use a gaussian deposition kernel of width rcl
ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL) # Delta function deposition kernel
# ds.eqsys.spi.setHeatAbsorbtion(SPI.HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS) # Remove all heat flowing through a disc 
#                                                                            of radius rcl from background plasma.
#                                                                            when assuming a local and immediate deposition
#                                                                            the heat absorbed by the ablated material is
#                                                                            immediately returned, so this term should then 
#                                                                            be switched off 
# ds.eqsys.spi.setCloudRadiusMode(SPI.CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT)

# Size of the neutral cloud (only relevant when using gaussian deposition or heat absorbtion)
rcl=0.01
ds.eqsys.spi.setRclPrescribedConstant(rcl)

# Set background ions
n_D = 1e20
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n_D)


# set collision settings
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
#ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_NEGLECT
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
#ds.collisions.lnlambda = Collisions.LNLAMBDA_CONSTANT
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL

# Kinetic grid settings
if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)
    ds.hottailgrid.setBiuniformGrid(psep=0.07,npsep=70)
    nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=radialgrid, T0=temperature.flatten())
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
    #ds.eqsys.f_hot.enableIonJacobian(False)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

######################
# Run the simulation #
######################
# Use the nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_LU)


ds.other.include('fluid', 'scalar')

filename_ending='nShardD'+str(nShardD)+'NinjD'+str(NinjD)+'nShardNe'+str(nShardNe)+'NinjNe'+str(NinjNe)+'vpD'+str(abs_vp_meanD)+'vpNe'+str(abs_vp_meanNe)+'hottail'+str(hotTailGrid_enabled)+'heat_transport'+str(use_heat_transport)+'f_hot_transport'+str(use_f_hot_transport)+'dBB'+str(dBOverB)	
folder_name='hottail_scan_cylindrical_DMS_profs/'

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)
	
if run_init:

	# Save settings to HDF5 file
	ds.save(folder_name+'init_settings'+filename_ending+'.h5')
	runiface(ds, folder_name+'output_init'+filename_ending+'.h5', quiet=False)


#######################
# RESTART set current #
#######################

do=DREAMOutput(folder_name+'output_init'+filename_ending+'.h5')
conductivity=do.other.fluid.conductivity.getData()
jprof=(1-(1-0.001**(1/0.41))*(radialgrid/radialgrid[-1])**2)**0.41
I0=15e6
j0=I0/np.trapz(2*np.pi*radialgrid*jprof,radialgrid)
efield=j0*jprof/conductivity[-1,:]

ds.eqsys.E_field.setPrescribedData(efield=efield, radius=radialgrid)

if run_init:
	# Save settings to HDF5 file
	ds.save(folder_name+'init_settings_'+filename_ending+'.h5')
	runiface(ds, folder_name+'output_init_'+filename_ending+'.h5', quiet=False)

##########################
# RESTART injection init #
##########################
# Used to accelerate the distribution to carry the right current

ds2 = DREAMSettings(ds)
ds2.fromOutput(folder_name+'output_init_'+filename_ending+'.h5')

ds2.timestep.setTmax(Tmax_init2)
ds2.timestep.setNt(Nt_init2)

if run_injection_init:
	ds2.save(folder_name+'init2_restart_settings_'+filename_ending+'.h5')
	runiface(ds2, folder_name+'output_init2_'+filename_ending+'.h5', quiet=False)



#####################
# RESTART injection #
#####################
# Here we actually make the first injection
ds3 = DREAMSettings(ds2)

# From now on, the temperature and electric field should be calculated self-consistently
if T_selfconsistent:
	ds3.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)
	ds3.eqsys.T_cold.setRecombinationRadiation(T_cold.RECOMBINATION_RADIATION_INCLUDED)
	
ds3.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds3.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall_R0 = E_wall*2*np.pi)

ds3.fromOutput(folder_name+'output_init2_'+filename_ending+'.h5',ignore=['v_p'])

# Now make the shards from the first injection start moving
ds3.eqsys.spi.setShardVelocitiesUniform(nShard=None,abs_vp_mean=abs_vp_meanD,abs_vp_diff=abs_vp_diffD,alpha_max=alpha_maxD,nDim=2,add=False, shards=slice(0,nShardD))

# If the first injection contains any impurities, the pressure and current will be 
# significantly perturbed already during this injection stage, and then we
# turn on the transport when the shards reach the plasma boundary
if use_heat_transport and molarFractionNe>0:
	print('Turn on transport already during first injection')
	dBB_mat=np.sqrt(R)*dBOverB*np.vstack((np.zeros((2,Nr)),np.ones((2,Nr))))
	t_edge=(radius_wall-radius[1])/np.max(-ds3.eqsys.spi.vp)
	print(t_edge)
	ds3.eqsys.T_cold.transport.setMagneticPerturbation(dBB=dBB_mat,r=radialgrid,t=[0,t_edge-1e-8,t_edge,1])
	ds3.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)
	
	if use_f_hot_transport and hotTailGrid_enabled:
		ds3.eqsys.n_re.transport.prescribeDiffusion(drr=sp.constants.pi*sp.constants.c*R*dBOverB**2*np.vstack((np.zeros((2,Nr)),np.ones((2,Nr)))),r=radialgrid,t=[0,0.00015999,0.00016,1])
		ds3.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
		ds3.eqsys.f_hot.transport.setMagneticPerturbation(dBB=np.sqrt(R)*dBOverB*np.ones(radialgrid.shape).reshape(1,-1),r=radialgrid,t=[0])
		ds3.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)

ds3.timestep.setTmax(Tmax_injection)
ds3.timestep.setNt(Nt_injection)
ds3.timestep.setNumberOfSaveSteps(int(Tmax_injection/1e-4))

if run_injection:
	ds3.save(folder_name+'injection_restart_settings_'+filename_ending+'.h5')
	runiface(ds3, folder_name+'output_restart_injection_'+filename_ending+'.h5', quiet=False)
	
##############
# Restart CQ #
##############
ds4 = DREAMSettings(ds3)
ds4.fromOutput(folder_name+'output_restart_injection_'+filename_ending+'.h5',ignore=['x_p','v_p'])
ds4.timestep.setTmax(Tmax_CQ)
ds4.timestep.setNt(Nt_CQ)
ds4.timestep.setNumberOfSaveSteps(int(Tmax_CQ/1e-4))

# Turn on transport, if any
if use_heat_transport:
	ds4.eqsys.T_cold.transport.setMagneticPerturbation(dBB=np.sqrt(R)*dBOverB*np.ones(radialgrid.shape).reshape(1,-1),r=radialgrid,t=[0])
	ds4.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)
if use_f_hot_transport and hotTailGrid_enabled:
	ds4.eqsys.n_re.transport.prescribeDiffusion(drr=sp.constants.pi*sp.constants.c*R*dBOverB**2*np.ones(radialgrid.shape).reshape(1,-1),r=radialgrid,t=[0])
	ds4.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
	ds4.eqsys.f_hot.transport.setMagneticPerturbation(dBB=np.sqrt(R)*dBOverB*np.ones(radialgrid.shape).reshape(1,-1),r=radialgrid,t=[0])
	ds4.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)


if run_CQ:
	ds4.save(folder_name+'CQ_restart_settings_'+filename_ending+'.h5')
	
	# Set the shards of the second injection into motion and advance them until
	# the fastest shards reach the plasma edge
	do4=DREAMOutput(folder_name+'output_restart_injection_'+filename_ending+'.h5')
	ds4.eqsys.spi.vp=do4.eqsys.v_p.data[-1,:].flatten()
	ds4.eqsys.spi.xp=do4.eqsys.x_p.data[-1,:].flatten()
	if(nShardNe>0):
		ds4.eqsys.spi.setShardVelocitiesUniform(nShard=None,abs_vp_mean=abs_vp_meanNe,abs_vp_diff=abs_vp_diffNe,alpha_max=alpha_maxNe,nDim=2,add=False, shards=slice(-nShardNe,None))

		t_edge=(radius_wall-radius[1])/np.max(-ds4.eqsys.spi.vp[-3*nShardNe::3])
		ds4.eqsys.spi.xp[-3*nShardNe:]=ds4.eqsys.spi.xp[-3*nShardNe:]+ds4.eqsys.spi.vp[-3*nShardNe:]*t_edge
	
	runiface(ds4, folder_name+'output_restart_CQ_'+filename_ending+'.h5', quiet=False)
	
#################
# Restart CQ 2+ #
#################
ds5 = DREAMSettings(ds4)
for iCQ in range(nCQ_restart_start-2,nCQ_restart):
	if iCQ==-1:
		ds5.fromOutput(folder_name+'output_restart_CQ_'+filename_ending+'.h5')
	else:
		ds5.fromOutput(folder_name+'output_restart_CQ'+str(iCQ+2)+'_'+filename_ending+'.h5')
	ds5.timestep.setTmax(Tmax_CQ_restart)
	ds5.timestep.setNt(Nt_CQ_restart)
	ds5.timestep.setNumberOfSaveSteps(int(Tmax_CQ_restart/1e-4))

	if iCQ==-1:
		do5=DREAMOutput(folder_name+'output_restart_CQ_'+filename_ending+'.h5')
	else:
		do5=DREAMOutput(folder_name+'output_restart_CQ'+str(iCQ+2)+'_'+filename_ending+'.h5')
	ds5.eqsys.spi.setInitialData(rp=do5.eqsys.Y_p.calcRadii(t=-1).flatten(),xp=do5.eqsys.x_p.data[-1,:].flatten(),vp=do5.eqsys.v_p.data[-1,:].flatten())
	ds5.save(folder_name+'CQ'+str(iCQ+3)+'_restart_settings_'+filename_ending+'.h5')
	runiface(ds5, folder_name+'output_restart_CQ'+str(iCQ+3)+'_'+filename_ending+'.h5', quiet=False)
	

