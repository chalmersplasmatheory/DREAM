#!/usr/bin/env python3
#
# This example validates the Dreicer runaway generation rate in DREAM.
# The simulation uses constant temperature, density and electric field,
# and generates runaway current through Dreicer mechanism only.
# Avalanche generation is disabled.
#
# Run as
#
#   $ ./generateD.py
#

import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.RadialGrid as RGrid
from DREAM import runiface

ds = DREAMSettings()

# Physical parameters (参考2kineticN示例)
E = 6       # Electric field strength (V/m) - 降低电场强度
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)

# Grid parameters
pMax = 1    # maximum momentum in units of m_e*c
Np   = 300  # 增加动量网格点数以提高精度
Nxi  = 20   # 增加pitch网格点数
tMax = 1e-3 # 延长模拟时间 (参考2kineticN)
Nt   = 200  # 增加时间步数

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions (fully ionized deuterium)
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Enable Dreicer generation using neural network model (most accurate)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Disable avalanche generation explicitly
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

# Disable Compton and tritium generation
ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_NEGLECT)
ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_NEGLECT)

# Collision settings (参考2kineticN)
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED 
ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS 

# Hot-tail grid settings
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)
ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)
ds.eqsys.f_hot.setParticleSource(particleSource=FHot.PARTICLE_SOURCE_IMPLICIT)

# Disable runaway grid (我们专注于流体逃逸电子密度)
ds.runawaygrid.setEnabled(False)

# Set up radial grid (参考2kineticN)
ds.radialgrid.setB0(5)           # Tesla
ds.radialgrid.setMinorRadius(0.22) # meters
ds.radialgrid.setWallRadius(0.22)  # meters
ds.radialgrid.setNr(1)           # Single radial point

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.preconditioner.setEnabled(False)
ds.solver.setVerbose(False)

# Include quantities needed for analysis
ds.other.include('fluid', 'nu_s', 'nu_D', 'lnLambda')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

# Output settings
ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('../outputs/dreicer_output.h5')

# Save settings to HDF5 file
ds.save('dreicer_settings.h5')

# Run the simulation
print("Settings saved to dreicer_settings.h5")
print("Running simulation...")
print(f"Parameters: E={E} V/m, n={n:.2e} m^-3, T={T} eV")
print(f"Grid: Np={Np}, Nxi={Nxi}, tMax={tMax}, Nt={Nt}")
do = runiface(ds, '../outputs/dreicer_output.h5', quiet=False)

print("Dreicer validation simulation completed.")
print("Output saved to ../outputs/dreicer_output.h5")
print("Settings saved to dreicer_settings.h5")