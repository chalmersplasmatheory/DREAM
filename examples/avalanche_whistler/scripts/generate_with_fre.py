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
#   $ ./generate_with_fre.py --a 0.3 --R 1.67 --E 0.05 --n 5e18 --T 2165
#

import argparse
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

# ── Command-line argument parsing ──────────────────────────────────────────────
parser = argparse.ArgumentParser(description='DREAM Dreicer + avalanche simulation with f_re grid.')
parser.add_argument('--E', type=float, default=0.05,
                    help='Electric field strength in V/m, default: 0.05')
parser.add_argument('--n', type=float, default=5e18,
                    help='Electron density in m^-3, default: 5e18')
parser.add_argument('--T', type=float, default=2165,
                    help='Temperature in eV, default: 2165')
parser.add_argument('--a', type=float, default=0.4,
                    help='Minor radius in meters, default: 0.4')
parser.add_argument('--R', type=float, default=1.67,
                    help='Major radius in meters, default: 1.67')
parser.add_argument('--B0', type=float, default=1.4,
                    help='On-axis magnetic field in Tesla, default: 1.4')
parser.add_argument('--Np-hot', type=int, default=100, dest='Np_hot',
                    help='Hot-tail momentum grid points, default: 100')
parser.add_argument('--Np-re', type=int, default=200, dest='Np_re',
                    help='Runaway momentum grid points, default: 200')
parser.add_argument('--Nxi', type=int, default=40,
                    help='Pitch grid points, default: 40')
parser.add_argument('--tMax', type=float, default=3.0,
                    help='Simulation time in seconds, default: 3.0')
parser.add_argument('--Nt', type=int, default=3000,
                    help='Number of time steps, default: 3000')
parser.add_argument('--output', type=str, default='../outputs/dreicer_with_fre_output.h5',
                    help='Output HDF5 file path, default: ../outputs/dreicer_with_fre_output.h5')
parser.add_argument('--source', type=str, default='off', choices=['on', 'off'],
                    help='Enable (kinetic) or disable avalanche source, default: off')
args = parser.parse_args()

ds = DREAMSettings()

# Physical parameters
E = args.E       # Electric field strength (V/m)
n = args.n       # Electron density (m^-3)
T = args.T       # Temperature (eV)

# Grid parameters
pMax_hot = 1     # maximum momentum in units of m_e*c for hot-tail grid
pMax_re = 50     # maximum momentum in units of m_e*c for runaway grid
Np_hot = args.Np_hot
Np_re  = args.Np_re
Nxi    = args.Nxi
tMax   = args.tMax
Nt     = args.Nt

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions (fully ionized deuterium)
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable Dreicer generation (kinetic simulation naturally captures Dreicer through distribution evolution)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

# Enable/disable avalanche generation (kinetic model) — controlled by --source
if args.source == 'on':
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=2.0)
else:
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
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)  # Enable synchrotron radiation
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)
ds.eqsys.f_hot.setParticleSource(particleSource=FHot.PARTICLE_SOURCE_IMPLICIT)

# Set initial runaway electron distribution to zero
ds.eqsys.f_re.setInitialProfiles(n0=0, T0=T)  # Initially no runaway electrons
ds.eqsys.f_re.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
ds.eqsys.f_re.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)  # Enable synchrotron radiation
ds.eqsys.f_re.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

# Set up radial grid
a  = args.a
R  = args.R
B0 = args.B0
ds.radialgrid.setType(RGrid.TYPE_ANALYTIC_TOROIDAL)
ds.radialgrid.setB0(B0)              # Tesla
ds.radialgrid.setShaping(psi=0.0001, GOverR0=B0)
ds.radialgrid.setMinorRadius(a)    # meters (small radius a)
ds.radialgrid.setMajorRadius(R)    # meters (major radius R0)
ds.radialgrid.setWallRadius(a)     # meters (wall radius = minor radius)
ds.radialgrid.setNr(1)                # Single radial point
nr=1
r_f = np.linspace(a*0.99, a, nr+1)
r_f *= a/r_f[-1]
ds.radialgrid.setCustomGridPoints(r_f)

# Set solver type
ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.preconditioner.setEnabled(True)  # Enable preconditioner for better convergence
ds.solver.setVerbose(True)  # Enable verbose output to see what's happening

# Include quantities needed for analysis
ds.other.include('fluid', 'nu_s', 'nu_D', 'lnLambda')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

# Output settings
ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename(args.output)

# Save settings to HDF5 file
ds.save('dreicer_with_fre_settings.h5')

# Run the simulation
print("Settings saved to dreicer_with_fre_settings.h5")
print("Running simulation...")
print(f"Parameters: E={E} V/m, n={n:.2e} m^-3, T={T} eV, a={a} m, R={R} m, B0={B0} T")
print(f"Grid: Np_hot={Np_hot}, Np_re={Np_re}, Nxi={Nxi}, tMax={tMax}, Nt={Nt}")
do = runiface(ds, args.output, quiet=False)

print("Dreicer validation simulation with f_re completed.")
print(f"Output saved to {args.output}")
print("Settings saved to dreicer_with_fre_settings.h5")