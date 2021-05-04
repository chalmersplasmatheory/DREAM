#!/usr/bin/env python3

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

np.random.seed(1)

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
run_init=False
run_injection_init=False
run_injection=True
run_CQ=True

nCQ_restart_start=2
nCQ_restart=0

use_heat_transport=True
use_f_hot_transport=False
transport_CQ_only=False
dBOverB=1e-3
#dBOverB=0

# time resolution of restarted simulation
Tmax_restart = 3.4e-3 # simulation time in seconds
Nt_restart = 340    # number of time steps
#Tmax_restart = 1e-3 # simulation time in seconds
#Nt_restart = 500    # number of time steps

Tmax_CQ = 8e-3
Nt_CQ = 1000
Tmax_CQ_restart = 0.2e-3
Nt_CQ_restart = 100


n_D = 1e20
#n_D = 5.3e19
n_D_inj_spi = 1e0

B0 = 5.3            # magnetic field strength in Tesla
E_initial = 0.00032 # initial electric field in V/m
E_wall = 0.0        # boundary electric field in V/m
# NOTE: it does not work to have self-consistent E-field with prescribed BC with E_wall=0, 
# since that leads to Psi_wall=0 constantly, which does not work when you have a relative tolerance
T_initial = 20e3    # initial temperature in eV
#T_initial = 5e3    # initial temperature in eV
#T_initial = 23e3    # initial temperature in eV

Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps

Tmax_init2 = 3e-3   # simulation time in seconds
Nt_init2 = 300         # number of time steps

Nr = 11             # number of radial grid points
Np = 120            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 3            # maximum momentum in m_e*c
times  = [0]        # times at which parameters are given
radius = [0, 2]     # span of the radial grid
dr=(radius[1]-radius[0])/(Nr+1)
radialgrid = np.linspace(radius[0]+dr/2,radius[-1]-dr/2,Nr)
radius_wall = 2.15  # location of the wall 

T_selfconsistent    = True
hotTailGrid_enabled = False

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)
ds.radialgrid.setWallRadius(radius_wall)

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)


# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
ds.eqsys.E_field.setBoundaryCondition()

# Set runaway generation rates
# ds.eqsys.n_re.setCompton(RE.COMPTON_RATE_ITER_DMS)
# ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)

temp_prof=(1-0.99*(radialgrid/radialgrid[-1])**2).reshape(1,-1)
#temp_prof=0.02+1*(radialgrid/radialgrid[-1])
#temp_prof=((1-0.9*(radialgrid/radialgrid[-1])**2)**2).reshape(1,-1)
temperature = T_initial*temp_prof
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radialgrid)
#temperature = T_initial * np.ones((len(times), len(radius)))
#ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)


#nShardD=1000
nShardD=1742
pelletDensity=205.9
pelletMolarMass=0.0020141
N_Avogadro=6.022e23


#Ninj=4.4e24
Ninj=2e24
kp=(Ninj/(6*np.pi**2*pelletDensity/pelletMolarMass*N_Avogadro*nShardD))**(-1/3)
kp=round(kp,4)
# kp=1137
def rp_distr(rp):
    return kn(0,rp*kp)*kp**2*rp

def sample_rp_distr(N):
    rp_integrate=np.linspace(1e-10,0.010,5000)
    cdf=integrate.cumtrapz(y=rp_distr(rp_integrate),x=rp_integrate)
    print(np.max(cdf))
    return np.interp(np.random.uniform(size=N),np.hstack((0.0,cdf)),rp_integrate)


#rp_init=0.002*np.ones(nShardD)
#rp_init=sample_rp_distr(nShardD)**(5/3)
rp_init=sample_rp_distr(nShardD)

#Ninj_obtained=np.sum(4*np.pi*rp_init**(9/5)/3*pelletDensity/pelletMolarMass*N_Avogadro)
#rp_init*=(Ninj/Ninj_obtained)**(5/9)
Ninj_obtained=np.sum(4*np.pi*rp_init**(3)/3*pelletDensity/pelletMolarMass*N_Avogadro)
rp_init*=(Ninj/Ninj_obtained)**(1/3)

print(np.sum(4*np.pi*rp_init**3/3*pelletDensity/pelletMolarMass*N_Avogadro))
print(kp)

# L0=0.1
xp_init=np.tile(np.array([radius_wall,0,0]),nShardD)
# xp_init=np.zeros(3*nShardD)
# xp_init[0::3]=radius[-1]+L0*np.random.uniform(size=nShardD)

print(xp_init[0])

alpha_max=0.17
abs_vp_mean=800
abs_vp_diff=0.2*abs_vp_mean
abs_vp_init=(abs_vp_mean+abs_vp_diff*(-1+2*np.random.uniform(size=nShardD)))
alpha=alpha_max*(-1+2*np.random.uniform(size=nShardD))
vp_init=np.zeros(3*nShardD)
vp_init[0::3]=-abs_vp_init*np.cos(alpha)
vp_init[1::3]=abs_vp_init*np.sin(alpha)

nShardNe=50
#nShardNe=0
pelletDensityNe=1444
pelletMolarMassNe=0.020183
N_Avogadro=6.022e23

if(nShardNe>0):
	NinjNe=10e24
	kpNe=(NinjNe/(6*np.pi**2*pelletDensityNe/pelletMolarMassNe*N_Avogadro*nShardNe))**(-1/3)
	kpNe=round(kpNe,4)
	alpha_maxNe=0.17
	abs_vp_meanNe=200
	abs_vp_diffNe=0.2*abs_vp_meanNe
	def rp_distr(rp):
		return kn(0,rp*kpNe)*kpNe**2*rp
    
if(nShardNe==1):
    rp_init_Ne=np.array([(NinjNe*3/(4*np.pi))**(1/3)])**(5/3)
    vp_initNe=np.array([-abs_vp_meanNe,0,0])
elif(nShardNe>1):
    rp_init_Ne=sample_rp_distr(nShardNe)**(5/3)
    NinjNe_obtained=np.sum(4*np.pi*rp_init_Ne**(9/5)/3*pelletDensityNe/pelletMolarMassNe*N_Avogadro)
    rp_init_Ne*=(NinjNe/NinjNe_obtained)**(5/9)
    abs_vp_initNe=(abs_vp_meanNe+abs_vp_diffNe*(-1+2*np.random.uniform(size=nShardNe)))
    alphaNe=alpha_maxNe*(-1+2*np.random.uniform(size=nShardNe))
    vp_initNe=np.zeros(3*nShardNe)
	#actually set the velocity of the neon shards when injected...
    
if(nShardNe>0):
	xp_init_Ne=np.tile(np.array([2.15,0,0]),nShardNe)
	rp_init=np.concatenate((rp_init,rp_init_Ne))
	vp_init=np.concatenate((vp_init, vp_initNe))
	xp_init=np.concatenate((xp_init,xp_init_Ne))
	nShard=nShardD+nShardNe

	SPIMolarFractionNe=np.zeros(nShard)
	SPIMolarFractionNe[nShardD:nShard]=1
	SPIMolarFractionD=np.ones(nShard)
	SPIMolarFractionD[nShardD:nShard]=1e-10
else:
	nShard=nShardD
	SPIMolarFractionNe=0.05*np.ones(nShard)
	SPIMolarFractionD=0.98*np.ones(nShard)


ds.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=np.zeros(3*nShard))

# Set ions
#n_D_prof=(1-0.9*(radialgrid/radialgrid[-1])**2)**(2/3)
density_D = n_D
density_D_inj = n_D_inj_spi

ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=density_D)
ds.eqsys.n_i.addIon(name='D_inj', Z=1, isotope=2, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=density_D_inj, SPIMolarFraction=SPIMolarFractionD)
ds.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e0, SPIMolarFraction=SPIMolarFractionNe)

R=6.2
kappa=1
ds.eqsys.spi.setVpVolNormFactor(R*kappa)

ds.eqsys.spi.setVelocity(SPI.VELOCITY_MODE_PRESCRIBED)
ds.eqsys.spi.setAblation(SPI.ABLATION_MODE_FLUID_NGS)
# ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL_GAUSSIAN)
ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL)
#ds.eqsys.spi.setHeatAbsorbtion(SPI.HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS)
#ds.eqsys.spi.setCloudRadiusMode(SPI.CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT)

rcl=0.01
ds.eqsys.spi.setRclPrescribedConstant(rcl)

if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)
    ds.hottailgrid.setBiuniformGrid(psep=0.07,npsep=70)
    nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    #ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=0, T0=T_initial)
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=radialgrid, T0=temperature.flatten())
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
    #ds.eqsys.f_hot.enableIonJacobian(False)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)


# Use the linear solver
# ds.solver.setType(Solver.LINEAR_IMPLICIT)

# Use the nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_LU)
ds.solver.setTolerance(reltol=0.001)
ds.solver.setMaxIterations(maxiter = 500)
#ds.solver.setVerbose(True)


ds.other.include('fluid', 'scalar')

#filename_ending='deposition'+str(ds.eqsys.spi.deposition)+'heatAbsorbtion'+str(ds.eqsys.spi.heatAbsorbtion)+'cloudRadiusMode'+str(ds.eqsys.spi.cloudRadiusMode)+'Nt'+str(Nt_restart)+'Nr'+str(Nr)
if(nShardNe>0):
	filename_ending='nShardD'+str(nShardD)+'kpD'+str(kp)+'nShardNe'+str(nShardNe)+'kpNe'+str(kpNe)+'vpD'+str(abs_vp_mean)+'vpNe'+str(abs_vp_meanNe)+'hottail'+str(hotTailGrid_enabled)+'heat_transport'+str(use_heat_transport)+'f_hot_transport'+str(use_f_hot_transport)+'dBB'+str(dBOverB)+'CQ_only'+str(transport_CQ_only)+'_new'
else:
	filename_ending='nShardD'+str(nShardD)+'kpD'+str(kp)+'vpD'+str(abs_vp_mean)+'_Ne_mixed'+str(SPIMolarFractionNe[0])+'_hottail'+str(hotTailGrid_enabled)+'heat_transport'+str(use_heat_transport)+'f_hot_transport'+str(use_f_hot_transport)+'dBB'+str(dBOverB)+'CQ_only'+str(transport_CQ_only)
	
folder_name='hottail_scan_cylindrical_DMS_profs/'
	
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
ds2.fromOutput(folder_name+'output_init_'+filename_ending+'.h5',ignore=['Y_p','x_p','v_p'])

#ds2.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
#ds2.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = E_wall*2*np.pi, wall_radius=radius_wall)


if T_selfconsistent:
	ds2.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

ds2.timestep.setTmax(Tmax_init2)
ds2.timestep.setNt(Nt_init2)
#ds2.timestep.setType(ttype=TimeStep.TYPE_ADAPTIVE)

if run_injection_init:
	ds2.save(folder_name+'init2_restart_settings_'+filename_ending+'.h5')
	runiface(ds2, folder_name+'output_init2_'+filename_ending+'.h5', quiet=False)



#####################
# RESTART injection #
#####################

ds3 = DREAMSettings(ds2)

ds3.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds3.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = E_wall*2*np.pi)

if use_heat_transport and (nShardNe<1 and not transport_CQ_only):
	dBB_mat=np.sqrt(R)*dBOverB*np.vstack((np.zeros((2,Nr)),np.ones((2,Nr))))
	#t_dBB=np.array([0,0.00015999,0.00016,1]).reshape(-1,1)*np.vstack((np.zeros((2,Nr)),np.ones((2,Nr))))
	ds3.eqsys.T_cold.transport.setMagneticPerturbation(dBB=dBB_mat,r=radialgrid,t=[0,0.00015999,0.00016,1])
	ds3.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)
	
	if use_f_hot_transport and hotTailGrid_enabled:
		ds3.eqsys.n_re.transport.prescribeDiffusion(drr=sp.constants.pi*sp.constants.c*R*dBOverB**2*np.vstack((np.zeros((2,Nr)),np.ones((2,Nr)))),r=radialgrid,t=[0,0.00015999,0.00016,1])
		ds3.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
		ds3.eqsys.f_hot.transport.setMagneticPerturbation(dBB=np.sqrt(R)*dBOverB*np.ones(radialgrid.shape).reshape(1,-1),r=radialgrid,t=[0])
		ds3.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)

# ds3.fromOutput(folder_name+'output_init2_'+filename_ending+'.h5',ignore=['Y_p','x_p','v_p'])
# ds3.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=vp_init)
ds3.fromOutput(folder_name+'output_init2_'+filename_ending+'.h5',ignore=['v_p'])
ds3.eqsys.spi.vp=vp_init

ds3.timestep.setTmax(Tmax_restart)
ds3.timestep.setNt(Nt_restart)
ds3.timestep.setNumberOfSaveSteps(int(Tmax_restart/1e-4))
#ds2.timestep.setType(ttype=TimeStep.TYPE_ADAPTIVE)

if run_injection:
	ds3.save(folder_name+'injection_restart_settings_'+filename_ending+'.h5')
	runiface(ds3, folder_name+'output_restart_injection_'+filename_ending+'.h5', quiet=False)
	
##############
# Restart CQ #
##############
ds4 = DREAMSettings(ds3)
ds4.fromOutput(folder_name+'output_restart_injection_'+filename_ending+'.h5',ignore=['Y_p','x_p','v_p'])
ds4.timestep.setTmax(Tmax_CQ)
ds4.timestep.setNt(Nt_CQ)
ds4.timestep.setNumberOfSaveSteps(int(Tmax_CQ/1e-4))

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
	
	do4=DREAMOutput(folder_name+'output_restart_injection_'+filename_ending+'.h5')
	vp=do4.eqsys.v_p.data[-1,:].flatten()
	xp=do4.eqsys.x_p.data[-1,:].flatten()
	if(nShardNe>0):
		vp_initNe[0::3]=-abs_vp_initNe*np.cos(alphaNe)
		vp_initNe[1::3]=abs_vp_initNe*np.sin(alphaNe)
		vp[3*nShardD:]=vp_initNe

		t_edge=0.15/np.max(-vp_initNe[0::3])
		xp_init_Ne[0::3]=2.15-t_edge*abs_vp_initNe*np.cos(alphaNe)
		xp_init_Ne[1::3]=t_edge*abs_vp_initNe*np.sin(alphaNe)
		xp[3*nShardD:]=xp_init_Ne
		
	rp=do4.eqsys.Y_p.calcRadii(t=-1).flatten()
	ds4.eqsys.spi.setInitialData(rp=rp,xp=xp,vp=vp)
	runiface(ds4, folder_name+'output_restart_CQ_'+filename_ending+'.h5', quiet=False)
	
#################
# Restart CQ 2+ #
#################
ds5 = DREAMSettings(ds4)
for iCQ in range(nCQ_restart_start-2,nCQ_restart):
	if iCQ==-1:
		ds5.fromOutput(folder_name+'output_restart_CQ_'+filename_ending+'.h5',ignore=['Y_p','x_p','v_p'])
	else:
		ds5.fromOutput(folder_name+'output_restart_CQ'+str(iCQ+2)+'_'+filename_ending+'.h5',ignore=['Y_p','x_p','v_p'])
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
	

