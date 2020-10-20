#!/usr/bin/env python3
#
# This script is intended to illustrate the energy balance by 
# plotting ohmic heating and radiative losses as a function of temperature
# at equilibrium ionization, similarly to figure 6 in Vallhagen et al JPP 2020. 
# This is achieved by setting a prescribed temperature profile at the values
# one wants to plot for and run a dynamic simulation until equilibrium has been reached
# (since equilibrium ionization settings does not seem to work yet).
#
# NOTE! Depending on the densities and temperatures one might have to adjust Tmax_restart_eq
# to be long enough to really reach sufficient equilibration!
#
# ###################################################################

import numpy as np
import sys
import matplotlib.pyplot as plt

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

# ds.eqsys.n_re.setEceff(Eceff=RE.COLLQTY_ECEFF_MODE_SIMPLE)

#############################
# Set simulation parameters #
#############################

n_D = 41e20 # deuterium density
n_Z = 0.08e20 # Impurity density

J=1.69e6 # Current density (For caculation of ohmic heating)

B0 = 5.3            # magnetic field strength in Tesla

Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps

Tmax_restart_ioniz = 2e-6
Nt_restart_ioniz = 500

Tmax_restart_eq = 30e-3
Nt_restart_eq = 1000

Tmax_restart_rad=1e-11
Nt_restart_rad=2

Nr = 151             # number of radial grid points
times  = [0]        # times at which parameters are given
radius = [0, 2]     # span of the radial grid
radialgrid = np.linspace(radius[0],radius[-1],Nr)
radius_wall = 2.15  # location of the wall 

E_initial = 0.001 # initial electric field in V/m (arbitrary value, does not affect the purpose of this script)
E_wall = 0.0001        # boundary electric field in V/m
# NOTE: it does not work to have self-consistent E-field with prescribed BC with E_wall=0, 
# since that leads to Psi_wall=0 constantly, which does not work when you have a relative tolerance

T_initial = np.logspace(np.log10(0.7),np.log10(2e3),Nr)    # initial temperature in eV

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)

# Set ions
Z0=1
Z=10

# If one wants to start from another initial ionization than fully ionized deuterium and neutral impurities
# Depending on the temperature range of interest, this can give a faster equilibration
"""
n_D_tmp=np.zeros(2)
n_D_tmp[0]=0*n_D
n_D_tmp[1]=1*n_D
n_D_tmp=n_D_tmp.reshape(-1,1)*np.ones((1,len(radius)))
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D_tmp,r=np.array(radius))

n_Z_tmp=np.zeros(Z+1)
n_Z_tmp[Z0]=n_Z
n_Z_tmp=n_Z_tmp.reshape(-1,1)*np.ones((1,len(radius)))
ds.eqsys.n_i.addIon(name='Ne', Z=Z, iontype=Ions.IONS_DYNAMIC, n=n_Z_tmp,r=np.array(radius))
"""

ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n_D)
ds.eqsys.n_i.addIon(name='Ne', Z=Z, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_Z)

# Since this script is intended to illustrate the energy balance at equilibrium ionization,
# it would be preferable to use these settings but that does not seem to work yet.
"""
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_EQUILIBRIUM, n=n_D)
ds.eqsys.n_i.addIon(name='Ne', Z=Z, iontype=Ions.IONS_EQUILIBRIUM, n=n_Z)
"""

temperature = T_initial * np.ones((len(times), len(radialgrid)))
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radialgrid)

# Set E_field 
efield = E_initial*np.ones((len(times), len(radius)))
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radius)
ds.eqsys.E_field.setBoundaryCondition(wall_radius=radius_wall)

# Disable runaway and hot-tail grid
ds.runawaygrid.setEnabled(False)
ds.hottailgrid.setEnabled(False)

# Use the nonlinear solver
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_LU)
ds.solver.setTolerance(reltol=0.01)
ds.solver.setMaxIterations(maxiter = 500)
#ds.solver.setVerbose(True)


ds.other.include('fluid', 'lnLambda','nu_s','nu_D')

# Save settings to HDF5 file
ds.save('init_settings.h5')
runiface(ds, 'output_init.h5', quiet=False)

#### Ionization #############
ds2 = DREAMSettings(ds)

ds2.timestep.setTmax(Tmax_restart_ioniz)
ds2.timestep.setNt(Nt_restart_ioniz)

ds2.save('ioniz_restart_settings.h5')

runiface(ds2, 'output_restart_ioniz.h5', quiet=False)

#### Equilibration ############
ds3 = DREAMSettings(ds2)

ds3.timestep.setTmax(Tmax_restart_eq)
ds3.timestep.setNt(Nt_restart_eq)

ds3.save('eq_restart_settings.h5')

runiface(ds3, 'output_restart_eq.h5', quiet=False)

#### Radiation ################
ds4 = DREAMSettings(ds3)

ds4.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

ds4.timestep.setTmax(Tmax_restart_rad)
ds4.timestep.setNt(Nt_restart_rad)

ds4.save('rad_restart_settings.h5')

runiface(ds4, 'output_restart_rad.h5', quiet=False)

################ Plot #################
do=DREAMOutput(ds4.output.filename)
sigma=do.other.fluid.conductivity[0,:]
rad=do.other.fluid.Tcold_radiation[0,:]
T=do.eqsys.T_cold[0,:]
plt.loglog(T,J**2/sigma/1e6)
plt.loglog(T,rad/1e6)
plt.show()
