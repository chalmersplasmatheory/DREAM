#!/usr/bin/env python3
#
# This example shows how to set up a pellet injection (by making an SPI with a single shard)
# in a JET-like setting, similar to the case studied in figure 1 in K Gál et al PPCF 2008
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
import DREAM.Settings.Equations.SPI as SPI


from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature

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

# time resolution of restarted simulation
Tmax_restart = 8e-3 # simulation time in seconds
Nt_restart = 200     # number of time steps

n_D = 2.8e19
n_D_inj_spi = 1e0

B0 = 5.3            # magnetic field strength in Tesla
E_initial = 0.00032 # initial electric field in V/m
E_wall = 0        # boundary electric field in V/m
T_initial=3.1e3

Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps
Nr = 50             # number of radial grid points
Np = 200            # number of momentum grid points
Nxi = 5             # number of pitch grid points
pMax = 1.0          # maximum momentum in m_e*c
times  = [0]        # times at which parameters are given
radius = [0, 1]     # span of the radial grid
radialgrid = np.linspace(radius[0],radius[-1],Nr)
# radius_wall = 2.15  # location of the wall 
radius_wall = 1.15  # location of the wall

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

# Set ions
n_D_prof=(1-0.9*(radialgrid/radialgrid[-1])**2)**(2/3)
density_D = n_D*n_D_prof
density_D_inj = n_D_inj_spi

ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=density_D, r=radialgrid)
ds.eqsys.n_i.addIon(name='D_inj', Z=1, isotope=2, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=density_D_inj, SPIMolarFraction=1.0)


# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
ds.eqsys.E_field.setBoundaryCondition()

# Set runaway generation rates
# ds.eqsys.n_re.setCompton(RE.COMPTON_RATE_ITER_DMS)
ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)

# Set temperature
temp_prof=((1-0.75*(radialgrid/radialgrid[-1])**2)**2).reshape(1,-1)
temperature = T_initial*temp_prof
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radialgrid)

# Pellet settings
rp_init=0.003
xp_init=[radius_wall,0,0]
vp_init=[-160,0,0]
R=2.96
ds.eqsys.spi.setInitialData(rp=rp_init,xp=xp_init,vp=[0,0,0]) # Set the velocity before the restart when the pellet should start to move
ds.eqsys.spi.setVpVolNormFactor(R)
ds.eqsys.spi.setVelocity(SPI.VELOCITY_MODE_PRESCRIBED)
ds.eqsys.spi.setAblation(SPI.ABLATION_MODE_FLUID_NGS)
ds.eqsys.spi.setDeposition(SPI.DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE)
ds.eqsys.spi.setHeatAbsorbtion(SPI.HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS)
ds.eqsys.spi.setCloudRadiusMode(SPI.CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT)

rcl=0.01
ds.eqsys.spi.setRclPrescribedConstant(rcl)

if not hotTailGrid_enabled:
    ds.hottailgrid.setEnabled(False)
else:
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)
    nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=0, T0=T_initial)
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)


# Use the nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_MKL)


ds.other.include('fluid', 'scalar')

filename_ending='deposition'+str(ds.eqsys.spi.deposition)+'heatAbsorbtion'+str(ds.eqsys.spi.heatAbsorbtion)+'cloudRadiusMode'+str(ds.eqsys.spi.cloudRadiusMode)+'Nr'+str(Nr)+'new'
folder_name='output/'

# Save settings to HDF5 file
ds.save(folder_name+'init_settings'+filename_ending+'.h5')
runiface(ds, folder_name+'output_init'+filename_ending+'.h5', quiet=False)


#######################
# RESTART set current #
#######################

do=DREAMOutput(folder_name+'output_init'+filename_ending+'.h5')
conductivity=do.other.fluid.conductivity.getData()
jprof=1
I0=2e6
j0=I0/np.trapz(2*np.pi*radialgrid*jprof,radialgrid)
efield=j0*jprof/conductivity[-1,:]

ds.eqsys.E_field.setPrescribedData(efield=efield, radius=radialgrid)

# Save settings to HDF5 file
ds.save(folder_name+'init_settings2'+filename_ending+'.h5')
runiface(ds, folder_name+'output_init2'+filename_ending+'.h5', quiet=False)


#####################
# RESTART injection #
#####################

ds2 = DREAMSettings(ds)

ds2.fromOutput(folder_name+'output_init2'+filename_ending+'.h5',ignore=['v_p'])
ds2.eqsys.spi.setInitialData(vp=vp_init)

ds2.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds2.eqsys.E_field.setBoundaryCondition(bctype = Efield.BC_TYPE_PRESCRIBED, inverse_wall_time = 0, V_loop_wall_R0 = E_wall*2*np.pi)

if T_selfconsistent:
    ds2.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

ds2.timestep.setTmax(Tmax_restart)
ds2.timestep.setNt(Nt_restart)
ds2.timestep.setNumberOfSaveSteps(int(Tmax_restart/1e-4))

ds2.save(folder_name+'injection_restart_settings'+filename_ending+'.h5')
runiface(ds2, folder_name+'output_restart_injection'+filename_ending+'.h5', quiet=False)


