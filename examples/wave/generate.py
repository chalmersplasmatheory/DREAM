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
E = .1   # Electric field strength (V/m)
n = 0.5e19    # Electron density (m^-3)
T = 1000     # Temperature (eV)

# Grid parameters
pMax = 15    # maximum momentum in units of m_e*c
Np   = 100   # number of momentum grid points
Nxi  = 30    # number of pitch grid points
tMax = 0.5   # simulation time in seconds
Nt   = 50   # number of time steps
Nr   = 3     # number of radial points

# Set up radial grid
ds.radialgrid.setB0(2.2)
ds.radialgrid.setMinorRadius(0.67)
ds.radialgrid.setWallRadius(1.00)
ds.radialgrid.setNr(Nr)

# Set E_field
t_E = np.array([0.0,0.45,0.46,0.5])
E = E*np.ones(4)
E[2:] = 1.
ds.eqsys.E_field.setPrescribedData(E, times=t_E)

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
    psep=1.0, npsep=int(Np/4),
    thetasep=np.pi/2,nthetasep_frac=.8)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)

# Set boundary condition type at pMax
#ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary

# wave settings (nt x nr)
r_wave = [0]
t_wave = [0.0, 0.2, 0.21, 0.39, 0.4, 0.5]
ppar_res = 3.5*np.ones((len(t_wave), len(r_wave))) # resonant momentum
Delta_ppar_res = 0.1*np.ones((len(t_wave), len(r_wave))) # width of resonant momentum
Dxx_int = np.zeros((len(t_wave), len(r_wave))) # strength of resonance
Dxx_int[2,:] = 0.01 # start of ramp
Dxx_int[3,:] = 0.01 # end of ramp
# Wave mode
ds.radialgrid.setWave(ppar_res, Delta_ppar_res, Dxx_int, r=r_wave, t=t_wave)
ds.eqsys.f_hot.setWaveMode(DistFunc.WAVE_MODE_GAUSSIAN)
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set solver type
# ds.solver.setType(Solver.NONLINEAR) # semi-implicit time stepping
ds.solver.setType(Solver.LINEAR_IMPLICIT) # time stepping
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

