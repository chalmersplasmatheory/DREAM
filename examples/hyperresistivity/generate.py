#!/usr/bin/env python3
#
# A simple example of the hyperrestive diffusion term in the Vloop equation.
#

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import sys

sys.path.append('../../py')
from DREAM import DREAMSettings
from DREAM import runiface

import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Solver as Solver


# Plasma parameters
ne = 5e19   # m^-3
Te = 4e3    # eV

# Tokamak parameters
a  = 0.5        # minor radius (m)
b  = a * 1.2    # wall radius (m)
B0 = 3.2        # on-axis magnetic field (m)
R0 = 1.65       # major radius (m)

# Initial current profile
j  = lambda r : (1 - (r/a)**4)**1.5   # target current density profile
Ip = 800e3      # target plasma current (A)

# Hyperresistive diffusion coefficient
mu0 = scipy.constants.mu_0
iNt = (2*np.pi*R0/a)**2 * 1e-3
mD  = scipy.constants.m_p + scipy.constants.m_n
VA = B0 / np.sqrt(mu0 * ne * mD)
#L0 = 1/144 * mu0 / (4*np.pi) * VA * iNt
L0 = 1e-4

# Innermost flux surface on which hyperresistive
# diffusion is to be enabled
rLmax   = 0.5*a

# Radial grid resolution
nr      = 40

# Width of region around rLmax which should be extra resolved
drLmax  = 0.1*a
# Number of grid points in the high resolution region
nrLmax  = 20

###############################################
# 1. PRIMARY SETTINGS (CALCULATE CONDUCTIVITY)

ds = DREAMSettings()

# Radial grid
ds.radialgrid.setType(RadialGrid.TYPE_CYLINDRICAL)
#ds.radialgrid.setNr(20)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(b)
ds.radialgrid.setB0(B0)

# Construct a custom radial grid with additional resolution
# just within rLmax...
if nrLmax >= nr-1:
    raise Exception('nrLmax must be strictly smaller than nr.')

rl  = rLmax-drLmax

nr1 = int(np.ceil((nr-nrLmax)*(rLmax/a)))
nr3 = int(np.ceil((nr-nrLmax)*(1-rLmax/a)))

r1  = np.linspace(0, rl, nr1+1)[:-1]
r2  = np.linspace(rl, rLmax, nrLmax)
r3  = np.linspace(rLmax, a, nr3+1)[1:]

r_f = np.concatenate((r1, r2, r3))
ds.radialgrid.setCustomGridPoints(r_f)

# Disable kinetic grids
ds.hottailgrid.setEnabled(False)
ds.runawaygrid.setEnabled(False)

# Constant plasma parameters
ds.eqsys.E_field.setPrescribedData(0)
ds.eqsys.T_cold.setPrescribedData(Te)
ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=ne)

ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
ds.solver.setType(Solver.NONLINEAR)

ds.timestep.setTmax(1e-6)
ds.timestep.setNt(1)

ds.other.include('fluid/conductivity')

# Calculate conductivity
do = runiface(ds, 'output_conductivity.h5', quiet=False)


###################################################
# 2. SET ELECTRIC FIELD (AND THUS CURRENT PROFILE)
jprof = j(do.grid.r)
j0 = Ip * 2*np.pi / do.grid.integrate(jprof)
E = j0*jprof / do.other.fluid.conductivity[-1,:]

ds.eqsys.E_field.setPrescribedData(E, radius=do.grid.r)

do = runiface(ds, 'output_init.h5', quiet=False)


#################################
# 3. APPLY HYPERRESISTIVE TERM
dsMain = DREAMSettings(ds)

tLambda = np.linspace(0, 2e-3, 10)
#rLambda = np.array([0, a])
rLambda = np.linspace(0, a, 20)

Lambda  = L0 * np.ones((tLambda.size, rLambda.size))
# Limit diffusion in time
Lambda[np.where(tLambda>1e-3),:] = 0
# Limit diffusion in space
Lambda[:,np.where(rLambda<rLmax)] = 0

dsMain.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
dsMain.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=0, R0=R0)

dsMain.eqsys.psi_p.setHyperresistivity(Lambda, radius=rLambda, times=tLambda)

dsMain.timestep.setTmax(3e-3)
dsMain.timestep.setNt(100)

dsMain.save('settings_main.h5')

do = runiface(dsMain, 'output_main.h5', quiet=False)

jmax = np.amax(do.eqsys.j_tot[-1,:])*1.1
do.eqsys.j_tot.plot(t=[0,-1], show=False)
plt.plot([rLmax, rLmax], [0, jmax], '--')
plt.plot([rl, rl], [0, jmax], 'k:')

plt.xlim([0, a])
plt.ylim([0, jmax])

plt.show()

