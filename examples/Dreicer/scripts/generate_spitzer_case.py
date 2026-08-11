#!/usr/bin/env python3
#
# This example simulates the electron distribution evolution under a constant 
# electric field for Spitzer conductivity calculation.
# Parameters:
# E = 0.4 V/m (E/Ec = 9.6 and E/ED = 0.056)
# T = 3 keV, n = 5e19 m^-3, Zeff = 2, B = 0

import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.RadialGrid as RGrid
from DREAM import runiface
from DREAM.Formulas.PlasmaParameters import getTauEETh, getEc, getED

ds = DREAMSettings()

# Physical parameters
E = 0.4     # Electric field strength (V/m)
T = 3e3     # Temperature (eV)
n = 5e19    # Electron density (m^-3)
Zeff = 2    # Effective charge

# Calculate critical fields
Ec = getEc(T, n)  # Connor-Hastie field
ED = getED(T, n)  # Dreicer field

print(f"E/Ec = {E/Ec:.2f}")
print(f"E/ED = {E/ED:.3f}")

# Calculate thermal collision time
tauEETh = getTauEETh(T, n)
print(f"tauEETh = {tauEETh:.3e} s")

# Grid parameters
pMax = 1    # maximum momentum in units of m_e*c
Np   = 300  # momentum grid points
Nxi  = 20   # pitch grid points

# Runaway grid parameters
pMax_re = 10   # Maximum momentum for runaway grid
Np_re   = 100  # Momentum grid points for runaway grid
Nxi_re  = 20   # Pitch grid points for runaway grid

# Calculate simulation time to reach 500 thermal collision times
tMax = 500 * tauEETh
Nt   = 200  # time steps

print(f"Simulation time: {tMax:.3e} s")

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions (fully ionized with effective charge Zeff)
ds.eqsys.n_i.addIon(name='Ion', Z=Zeff, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable runaway generation (focus on Spitzer conductivity)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_NEGLECT)
ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_NEGLECT)
ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_NEGLECT)

# Collision settings
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
ds.collisions.lnlambda            = Collisions.LNLAMBDA_THERMAL

# Hot-tail grid settings
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

# Enable runaway grid
ds.runawaygrid.setEnabled(True)
ds.runawaygrid.setNxi(Nxi_re)
ds.runawaygrid.setNp(Np_re)
ds.runawaygrid.setPmax(pMax_re)

# Set up radial grid
ds.radialgrid.setB0(1e-6)        # Very small magnetic field (effectively B = 0)
ds.radialgrid.setMinorRadius(0.1) # Small radius
ds.radialgrid.setWallRadius(0.1) 
ds.radialgrid.setNr(1)           # Single radial point

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.solver.preconditioner.setEnabled(False)
ds.solver.setVerbose(False)

# Include quantities needed for analysis (remove tauEETh as it's not a valid other quantity)
ds.other.include('fluid', 'nu_s', 'nu_D', 'lnLambda')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

# Output settings
ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('../outputs/spitzer_output.h5')

# Save settings to HDF5 file
ds.save('spitzer_settings.h5')

# Run the simulation
print("Settings saved to spitzer_settings.h5")
print("Running simulation...")
print(f"Parameters: E={E} V/m, n={n:.2e} m^-3, T={T} eV, Zeff={Zeff}")
print(f"Grid: Np={Np}, Nxi={Nxi}, tMax={tMax:.2e}, Nt={Nt}")
do = runiface(ds, '../outputs/spitzer_output.h5', quiet=False)

print("Spitzer conductivity simulation completed.")
print("Output saved to ../outputs/spitzer_output.h5")
print("Settings saved to spitzer_settings.h5")