#!/usr/bin/env python3
#
# This example shows how to set up a simple test of the 
# implementation of prescribed ion diffusion coefficients in DREAM
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
from DREAM import runiface
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Solver as Solver

ds = DREAMSettings()

tMax_init = 1e-6 # simulation time in seconds
Nt_init   = 500   # number of time steps

tMax = 10e-3 # simulation time in seconds
Nt   = 200   # number of time steps

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)
ds.hottailgrid.setEnabled(False)

# Set up radial grid
Nr = 10
radius = [0, 2]     # span of the radial grid
ds.radialgrid.setB0(5.3)
ds.radialgrid.setMinorRadius(radius[-1])
ds.radialgrid.setWallRadius(2.15)
ds.radialgrid.setNr(Nr)

times  = [0]        # times at which parameters are given
dr=(radius[1]-radius[0])/(Nr+1)
radialgrid = np.linspace(radius[0]+dr/2,radius[-1]-dr/2,Nr)

# Physical parameters
E = 0       # Electric field strength (V/m)
n = 1e20 * np.ones(Nr)
#n = 1e20  * (1-0.99*radialgrid/radius[-1]) # Ion density (m^-3)s
n[5]=n[5]+5e19
charged_prescribed_diffusion = 10
neutral_prescribed_diffusion = 50
T = 1     # Temperature (eV)

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, 
    charged_diffusion_mode = Ions.ION_CHARGED_DIFFUSION_MODE_PRESCRIBED, charged_prescribed_diffusion =  charged_prescribed_diffusion,
    neutral_diffusion_mode = Ions.ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED, neutral_prescribed_diffusion =  neutral_prescribed_diffusion,
    n=n, r=radialgrid)

# Set solver type
#ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.setType(Solver.NONLINEAR)
ds.solver.preconditioner.setEnabled(False)

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax_init)
ds.timestep.setNt(Nt_init)

ds.output.setTiming(stdout=True, file=True)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

runiface(ds, 'output.h5', quiet=False)

ds2 = DREAMSettings(ds)
ds2.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

# Set time stepper
ds2.timestep.setTmax(tMax)
ds2.timestep.setNt(Nt)

runiface(ds2, 'output2.h5', quiet=False)



