#!/usr/bin/env python3
#
# This example validates the Dreicer runaway generation rate in DREAM.
# The simulation uses constant temperature, density and electric field,
# and generates runaway current through Dreicer mechanism only.
# Avalanche generation is disabled.
# This version enables the runaway electron grid to capture f_re data.
#
# Run as
#
#   $ ./generate_with_fre.py
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
E = 0.05    # Electric field strength (V/m)
n = 5e18    # Electron density (m^-3)
T = 2165    # Temperature (eV) - 2.165 keV

# Grid parameters
pMax_hot = 1    # maximum momentum in units of m_e*c for hot-tail grid
pMax_re = 10    # maximum momentum in units of m_e*c for runaway grid
Np_hot   = 100  # 动量网格点数（减少以加快计算）
Np_re   = 100   # 动量网格点数（减少以加快计算）
Nxi  = 40       # pitch网格点数（减少以加快计算）
tMax = 1e-1     # 模拟时间 (100 ms)
Nt   = 5000      # 时间步数

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions (fully ionized deuterium)
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Enable Dreicer generation using neural network model (most accurate)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Enable avalanche generation (kinetic model - most accurate)
# pCutAvalanche: momentum cutoff for avalanche generation (in units of m_e*c)
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=2.0)

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
ds.hottailgrid.setNp(Np_hot)
ds.hottailgrid.setPmax(pMax_hot)
ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)

# Runaway electron grid settings
ds.runawaygrid.setEnabled(True)
ds.runawaygrid.setNxi(Nxi)
ds.runawaygrid.setNp(Np_re)
ds.runawaygrid.setPmax(pMax_re)
ds.runawaygrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)
ds.eqsys.f_hot.setParticleSource(particleSource=FHot.PARTICLE_SOURCE_IMPLICIT)

# Set initial runaway electron distribution to zero
ds.eqsys.f_re.setInitialProfiles(n0=0, T0=T)  # Initially no runaway electrons
ds.eqsys.f_re.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
ds.eqsys.f_re.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)
ds.eqsys.f_re.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

# Set up radial grid
ds.radialgrid.setB0(1.4)              # Tesla
ds.radialgrid.setMinorRadius(0.67)    # meters (small radius a)
ds.radialgrid.setMajorRadius(1.67)    # meters (major radius R0)
ds.radialgrid.setWallRadius(0.67)     # meters (wall radius = minor radius)
ds.radialgrid.setNr(1)                # Single radial point

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
ds.output.setFilename('../outputs/dreicer_with_fre_output.h5')

# Save settings to HDF5 file
ds.save('dreicer_with_fre_settings.h5')

# Run the simulation
print("Settings saved to dreicer_with_fre_settings.h5")
print("Running simulation...")
print(f"Parameters: E={E} V/m, n={n:.2e} m^-3, T={T} eV")
print(f"Grid: Np_hot={Np_hot}, Np_re={Np_re}, Nxi={Nxi}, tMax={tMax}, Nt={Nt}")
do = runiface(ds, '../outputs/dreicer_with_fre_output.h5', quiet=False)

print("Dreicer validation simulation with f_re completed.")
print("Output saved to ../outputs/dreicer_with_fre_output.h5")
print("Settings saved to dreicer_with_fre_settings.h5")