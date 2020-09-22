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

n_D = 1e20
n_Z = 1e20

B0 = 5.3            # magnetic field strength in Tesla
T_initial = 5    # initial temperature in eV

Tmax_init = 1e-11   # simulation time in seconds
Nt_init = 2         # number of time steps
Nr = 151             # number of radial grid points
times  = [0]        # times at which parameters are given
radius = [0, 2]     # span of the radial grid
radialgrid = np.linspace(radius[0],radius[-1],Nr)
radius_wall = 2.15  # location of the wall 

E_initial = np.linspace(0,100,Nr) # initial electric field in V/m
E_wall = 0.0001        # boundary electric field in V/m
# NOTE: it does not work to have self-consistent E-field with prescribed BC with E_wall=0, 
# since that leads to Psi_wall=0 constantly, which does not work when you have a relative tolerance

# Set up radial grid
ds.radialgrid.setB0(B0)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setNr(Nr)

# Set time stepper
ds.timestep.setTmax(Tmax_init)
ds.timestep.setNt(Nt_init)

# Set ions
Z0=1
Z=18

n_D_tmp=np.zeros(2)
n_D_tmp[0]=0*n_D
n_D_tmp[1]=1*n_D
n_D_tmp=n_D_tmp.reshape(-1,1)*np.ones((1,len(radius)))
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D_tmp,r=np.array(radius))

n_Z_tmp=np.zeros(Z+1)
n_Z_tmp[Z0]=n_Z
n_Z_tmp=n_Z_tmp.reshape(-1,1)*np.ones((1,len(radius)))
ds.eqsys.n_i.addIon(name='Ar', Z=Z, iontype=Ions.IONS_DYNAMIC, n=n_Z_tmp,r=np.array(radius))

# If one wants to use equilibrium ionization (does not work yet)
"""
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_EQUILIBRIUM, n=n_D)
ds.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_EQUILIBRIUM, n=n_Z)
"""


# Set E_field 
efield = E_initial*np.ones((len(times), len(radialgrid)))
print(efield.shape)
ds.eqsys.E_field.setPrescribedData(efield=efield, times=times, radius=radialgrid)
ds.eqsys.E_field.setBoundaryCondition(wall_radius=radius_wall)

# Set runaway generation rates
ds.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID)

temperature = T_initial * np.ones((len(times), len(radius)))
ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times, radius=radius)

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

ds.output.setFilename('output_init.h5')

# Save settings to HDF5 file
ds.save('init_settings.h5')
runiface(ds, ds.output.filename, quiet=False)

ds2=DREAMSettings(ds)
ds2.eqsys.n_re.setAvalanche(RE.AVALANCHE_MODE_FLUID_HESSLOW)
runiface(ds2, 'output_init2.h5', quiet=False)

do=DREAMOutput('output_init.h5')
do2=DREAMOutput('output_init2.h5')
E=do.eqsys.E_field[0,:]
plt.plot(E,do.other.fluid.GammaAva[0,:]*1e-3,label='DREAM')
plt.plot(E,do2.other.fluid.GammaAva[0,:]*1e-3,label='GO-like')
plt.legend()
plt.xlabel('$E$ [V/m]')
plt.ylabel('$\Gamma$ [1/ms]')

plt.show()

