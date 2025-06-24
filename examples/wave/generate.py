#!/usr/bin/env python3
#
# This script generates a simple CODE-like Dreicer runaway simulation
# with additional pitch scattering due to wave-injection.
#
# Run as
#
#   $ ./basic.py
#   $ ../../build/iface/dreami dream_settings.h5
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

ds = DREAMSettings()

# Physical parameters
E = .2       # Electric field strength (V/m)
n = 1.4e19    # Electron density (m^-3)
T = 800     # Temperature (eV)

# Grid parameters
pMax = 50    # maximum momentum in units of m_e*c
Np   = 200   # number of momentum grid points
Nxi  = 30    # number of pitch grid points
tMax = 0.7   # simulation time in seconds
Nt   = 20    # number of time steps

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

# Hot-tail grid settings
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)
ds.hottailgrid.setBiuniformGrid(
    psep=0.5, npsep=50,
    thetasep=np.pi/2,nthetasep_frac=.8)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)

# Set boundary condition type at pMax
#ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary

# wave settings (nt x nr)
r_wave = [0]
t_wave = [0.0, 0.3, 0.4, 0.7]
ppar_res = 0.5*np.ones((len(t_wave), len(r_wave))) # resonant momentum
Delta_ppar_res = 0.05*np.ones((len(t_wave), len(r_wave))) # width of resonant momentum
Dxx_int = np.zeros((len(t_wave), len(r_wave))) # strength of resonance
Dxx_int[1:3,:] = 1E-3 # only turn on wave in certain time interval
# Wave mode
ds.radialgrid.setWave(ppar_res, Delta_ppar_res, Dxx_int, r=r_wave, t=t_wave)
ds.eqsys.f_hot.setWaveMode(DistFunc.WAVE_MODE_GAUSSIAN)
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
Nr = 3 # number of radial points
ds.radialgrid.setB0(1.4)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setWallRadius(0.22)
ds.radialgrid.setNr(Nr)

# Set solver type
ds.solver.setType(Solver.NONLINEAR) # semi-implicit time stepping
ds.solver.tolerance.set(reltol=1e-4)
ds.solver.setVerbose(False)

# include otherquantities to save to output
#ds.other.include('fluid/runawayRate','ripple')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('outputs.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

from DREAM import runiface
ds.solver.setVerbose(True)
do = runiface(ds)

