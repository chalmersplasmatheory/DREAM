#!/usr/bin/env python3
#
# This example shows how to set up a self-consistent fluid DREAM run,
# where no kinetic equations are solved, but the electric field and
# temperature are evolved self-consistently.
#
# Run as
#
#   $ ./basic.py
#
# ###################################################################

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


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

np.random.seed(1)

ds = DREAMSettings()

# set collision settings
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
#ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_NEGLECT
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
#ds.collisions.lnlambda = Collisions.LNLAMBDA_CONSTANT
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL

#############################
# Set simulation parameters #
#############################

# time resolution of restarted simulation
Tmax_restart = 26e-3 # simulation time in seconds
Nt_restart = 5000    # number of time steps

n_D = 1e20
#n_D = 5.3e19
n_D_inj_spi = 1e0

B0 = 5.3            # magnetic field strength in Tesla
E_initial = 0.00032 # initial electric field in V/m
E_wall = 0.0        # boundary electric field in V/m
# NOTE: it does not work to have self-consistent E-field with prescribed BC with E_wall=0, 
# since that leads to Psi_wall=0 constantly, which does not work when you have a relative tolerance
#T_initial = 20e3    # initial temperature in eV
#T_initial = 5e3    # initial temperature in eV
T_initial = 23e3    # initial temperature in eV

Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps

Tmax_init2 = 3e-4   # simulation time in seconds
Nt_init2 = 300         # number of time steps

Nr = 11             # number of radial grid points
Np = 60            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 3            # maximum momentum in m_e*c
times  = [0]        # times at which parameters are given
radius = [0, 2]     # span of the radial grid
radialgrid = np.linspace(radius[0],radius[-1],Nr)
radius_wall = 2.15  # location of the wall 

T_selfconsistent    = True
hotTailGrid_enabled = True

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)

# Set ions
#n_D_prof=(1-0.9*(radialgrid/radialgrid[-1])**2)**(2/3)
density_D = n_D
density_D_inj = n_D_inj_spi

ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=density_D)
ds.eqsys.n_i.addIon(name='D_inj', Z=1, isotope=2, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=density_D_inj, SPIMolarFraction=1.0)
#ds.eqsys.n_i.addIon(name='Ar', Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e10)


# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
ds.eqsys.E_field.setBoundaryCondition(wall_radius=radius_wall)

# Set runaway generation rates
# ds.eqsys.n_re.setCompton(RE.COMPTON_RATE_ITER_DMS)
# ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)

#temp_prof=(1-0.99*(radialgrid/radialgrid[-1])**2).reshape(1,-1)
#temp_prof=0.02+1*(radialgrid/radialgrid[-1])
temp_prof=((1-0.9*(radialgrid/radialgrid[-1])**2)**2).reshape(1,-1)
temperature = T_initial*temp_prof
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radialgrid)
#temperature = T_initial * np.ones((len(times), len(radius)))
#ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)


#nShard=1000
nShard=564
pelletDensity=205.9
pelletMolarMass=0.0020141
N_Avogadro=6.022e23


#Ninj=2.2e24
Ninj=5e24
kp=(Ninj/(6*np.pi**2*pelletDensity/pelletMolarMass*N_Avogadro*nShard))**(-1/3)
kp=round(kp,4)
# kp=1137
def rp_distr(rp):
    return kn(0,rp*kp)*kp**2*rp

def sample_rp_distr(N):
    rp_integrate=np.linspace(1e-10,0.010,5000)
    cdf=integrate.cumtrapz(y=rp_distr(rp_integrate),x=rp_integrate)
    print(np.max(cdf))
    return np.interp(np.random.uniform(size=N),np.hstack((0.0,cdf)),rp_integrate)


#rp_init=0.002*np.ones(nShard)
#rp_init=4*np.pi*sample_rp_distr(nShard)**3/3*pelletDensity/pelletMolarMass*N_Avogadro
rp_init=sample_rp_distr(nShard)**(5/3)

print(np.sum(4*np.pi*rp_init**3/3*pelletDensity/pelletMolarMass*N_Avogadro))
print(kp)

# L0=0.1
xp_init=np.tile(np.array([radius_wall,0,0]),nShard)
# xp_init=np.zeros(3*nShard)
# xp_init[0::3]=radius[-1]+L0*np.random.uniform(size=nShard)

print(xp_init[0])

alpha_max=0.17
abs_vp_mean=200
abs_vp_diff=0.2*abs_vp_mean
abs_vp_init=(abs_vp_mean+abs_vp_diff*(-1+2*np.random.uniform(size=nShard)))
alpha=alpha_max*(-1+2*np.random.uniform(size=nShard))
vp_init=np.zeros(3*nShard)
vp_init[0::3]=-abs_vp_init*np.cos(alpha)
vp_init[1::3]=abs_vp_init*np.sin(alpha)

ds.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=np.zeros(3*nShard))

R=6.2*1.8
ds.eqsys.spi.setVpVolNormFactor(R)

ds.eqsys.spi.setVelocity(SPI.VELOCITY_MODE_PRESCRIBED)
ds.eqsys.spi.setAblation(SPI.ABLATION_MODE_KINETIC_NGS)
# ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL_GAUSSIAN)
ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL)
# ds.eqsys.spi.setHeatAbsorbtion(SPI.HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS)
# ds.eqsys.spi.setCloudRadiusMode(SPI.CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT)

rcl=0.01
ds.eqsys.spi.setRclPrescribedConstant(rcl)

if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)
    #ds.hottailgrid.setBiuniformGrid(psep=0.06,npsep=50)
    nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    #ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=0, T0=T_initial)
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=0.99*nfree_initial, rT0=radialgrid, T0=temperature.flatten())
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


ds.other.include('fluid', 'lnLambda','nu_s','nu_D')

filename_ending='deposition'+str(ds.eqsys.spi.deposition)+'heatAbsorbtion'+str(ds.eqsys.spi.heatAbsorbtion)+'cloudRadiusMode'+str(ds.eqsys.spi.cloudRadiusMode)+'Nt'+str(Nt_restart)+'Nr'+str(Nr)+'superthermal'

# Save settings to HDF5 file
ds.save('init_settings'+filename_ending+'.h5')
runiface(ds, 'output_init'+filename_ending+'.h5', quiet=False)


#######################
# RESTART set current #
#######################

do=DREAMOutput('output_init'+filename_ending+'.h5')
conductivity=do.other.fluid.conductivity.getData()
jprof=(1-(1-0.001**(1/0.41))*(radialgrid/radialgrid[-1])**2)**0.41
# efield=1.81e6*jprof/conductivity[-1,:]
efield=1.69e6*jprof/conductivity[-1,:]

ds.eqsys.E_field.setPrescribedData(efield=efield, radius=radialgrid)

# Save settings to HDF5 file
ds.save('init_settings_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5')
runiface(ds, 'output_init_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5', quiet=False)

##########################
# RESTART injection init #
##########################

ds2 = DREAMSettings(ds)

ds2.fromOutput('output_init_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5',ignore=['r_p','x_p','v_p'])
ds.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=vp_init)

#ds2.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
#ds2.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = E_wall*2*np.pi, wall_radius=radius_wall)


if T_selfconsistent:
    ds2.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

ds2.timestep.setTmax(Tmax_init2)
ds2.timestep.setNt(Nt_init2)
#ds2.timestep.setType(ttype=TimeStep.TYPE_ADAPTIVE)

ds2.save('init2_restart_settings_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5')
runiface(ds2, 'output_init2_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5', quiet=False)



#####################
# RESTART injection #
#####################

ds3 = DREAMSettings(ds2)
ds3.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=vp_init)

ds3.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds3.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall = E_wall*2*np.pi, wall_radius=radius_wall)

ds3.fromOutput('output_init2_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5',ignore=['r_p','x_p','v_p'])
ds3.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=vp_init)

ds3.timestep.setTmax(Tmax_restart)
ds3.timestep.setNt(Nt_restart)
#ds2.timestep.setType(ttype=TimeStep.TYPE_ADAPTIVE)

ds3.save('injection_restart_settings_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5')
runiface(ds3, 'output_restart_injection_nShard'+str(nShard)+'kp'+str(kp)+''+filename_ending+'.h5', quiet=False)


