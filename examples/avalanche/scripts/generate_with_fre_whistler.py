#!/usr/bin/env python3
#
# This example tests quasilinear diffusion from whistler waves in DREAM.
# Based on the avalanche/generate_with_fre.py example, but adds wave-particle
# interactions using QUADRE parameters (476 MHz whistler wave).
#
# Run as
#
#   $ ./generate_with_fre_whistler.py [--amplitude AMPLITUDE] [--a A] [--R R]
#

import numpy as np
import sys
import argparse

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

# Parse command line arguments
parser = argparse.ArgumentParser(description='DREAM Quasilinear Diffusion Test - Whistler Wave')
parser.add_argument('--amplitude', type=float, default=1e3,
                    help='Wave amplitude (normalized units), default: 1e3')
parser.add_argument('--a', type=float, default=0.3,
                    help='Minor radius in meters, default: 0.3')
parser.add_argument('--R', type=float, default=1.67,
                    help='Major radius in meters, default: 1.67')
args = parser.parse_args()

ds = DREAMSettings()

# ============================================================================
# Physical parameters (matching QUADRE simulation)
# ============================================================================
E = 0.05 # Electric field strength (V/m)
n = 5e18           # Electron density (m^-3)
T = 2165           # Temperature (eV) - 2.165 keV
B0 = 1.4           # Magnetic field (Tesla)

# Grid parameters
pMax_hot = 1       # Maximum momentum for hot-tail grid
pMax_re = 50       # Maximum momentum for runaway grid (must be > pMax_hot)
Np_hot   = 100     # Momentum grid points
Np_re    = 200     # Runaway momentum grid points
Nxi      = 40      # Pitch-angle grid points
tMax     = 3.0     # Simulation time (seconds)
Nt       = 3000    # Time steps

print("="*70)
print("DREAM Quasilinear Diffusion Test - 476 MHz Whistler Wave")
print("="*70)
print(f"\nPlasma Parameters:")
print(f"  B0 = {B0} T")
print(f"  n_e = {n:.2e} m^-3")
print(f"  T = {T} eV")
print(f"  E = {E:.2f} V/m (E/Ec = 0.045)")

# ============================================================================
# Set electric field and temperature
# ============================================================================
ds.eqsys.E_field.setPrescribedData(E)
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions (fully ionized deuterium)
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# ============================================================================
# Runaway electron settings
# ============================================================================
# Disable Dreicer generation (kinetic simulation captures it naturally)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

# Enable avalanche generation (kinetic model)
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=2.0)

# Disable Compton and tritium generation
ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_NEGLECT)
ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_NEGLECT)

# Collision settings
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED 
ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS 

# ============================================================================
# Hot-tail grid settings
# ============================================================================
ds.hottailgrid.setNxi(Nxi)
ds.hottailgrid.setNp(Np_hot)
ds.hottailgrid.setPmax(pMax_hot)
ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)

# ============================================================================
# Runaway electron grid settings
# ============================================================================
ds.runawaygrid.setEnabled(True)
ds.runawaygrid.setNxi(Nxi)
ds.runawaygrid.setNp(Np_re)
ds.runawaygrid.setPmax(pMax_re)
ds.runawaygrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 0.1)

# ============================================================================
# Initial distribution and physics terms
# ============================================================================
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0)
ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)
ds.eqsys.f_hot.setParticleSource(particleSource=FHot.PARTICLE_SOURCE_IMPLICIT)

# Set initial runaway electron distribution to zero
ds.eqsys.f_re.setInitialProfiles(n0=0, T0=T)
ds.eqsys.f_re.setBoundaryCondition(DistFunc.BC_F_0)
ds.eqsys.f_re.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)
ds.eqsys.f_re.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

# ============================================================================
# QUASILINEAR DIFFUSION FROM WHISTLER WAVES (QUADRE parameters)
# ============================================================================
print(f"\nQuasilinear Diffusion Configuration:")
print(f"  Wave frequency: 476 MHz")
print(f"  k_parallel: -41 m^-1")

# QUADRE wave parameters from simulation.log
quadre_wave_params = {
    'k_main': 54.58,              # Total wavenumber (m^-1)
    'ktheta_main': 2.42,          # Propagation angle (rad) - obtuse angle for backward propagation
    'k_range': [52.61, 56.55],    # Wavenumber range (m^-1)
    'ktheta_range': [2.4, 2.45]  # Angle range (rad) - centered around 2.42 rad
}

print(f"  k range: [{quadre_wave_params['k_range'][0]:.2f}, {quadre_wave_params['k_range'][1]:.2f}] m^-1")
print(f"  θ_k range: [{quadre_wave_params['ktheta_range'][0]:.2f}, {quadre_wave_params['ktheta_range'][1]:.2f}] rad")

# Enable quasilinear diffusion on f_re (runaway electrons)
# This is the correct target because:
# - Whistler waves resonate with high-energy electrons (p ~ 16, ~7.7 MeV)
# - These are in the runaway energy range
# - Goal is to suppress runaway growth via wave-particle interactions
ds.eqsys.f_re.setQuasilinearDiffusion(
    enabled=True,
    quadre_params=quadre_wave_params,
    num_k=8,              # Match QUADRE discretization
    num_ktheta=20,        # Match QUADRE discretization
    amplitude=args.amplitude,      # Wave amplitude (normalized units) - PHYSICAL VALUE from QUADRE (δB/B₀ ~ 10⁻⁵)
    harmonic_mode='both',  # Include n = -2,-1,0,+1,+2
    use_simple_dispersion=True  # Use simplified dispersion for testing (ω = k|k_∥| * w)
)

print(f"  ✓ Quasilinear diffusion enabled")
print(f"  - Spectrum: uniform grid")
print(f"  - Modes: {8} × {20} = {8*20}")
print(f"  - Harmonics: n ∈ {{-2, -1, 0, +1, +2}}")

# ============================================================================
# Radial grid setup (copied from generate_with_fre.py)
# ============================================================================
a = args.a   # Minor radius (meters)
R = args.R   # Major radius (meters)

ds.radialgrid.setType(RGrid.TYPE_ANALYTIC_TOROIDAL)
ds.radialgrid.setB0(B0)
ds.radialgrid.setShaping(psi=0.0001, GOverR0=B0)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setMajorRadius(R)
ds.radialgrid.setWallRadius(a)
ds.radialgrid.setNr(1)

nr = 1
r_f = np.linspace(a*0.99, a, nr+1)   
r_f *= a/r_f[-1]
ds.radialgrid.setCustomGridPoints(r_f)

# ============================================================================
# Solver settings
# ============================================================================
ds.solver.setType(Solver.LINEAR_IMPLICIT)
ds.solver.preconditioner.setEnabled(True)
ds.solver.setVerbose(True)

# Include quantities for analysis
ds.other.include('fluid', 'nu_s', 'nu_D', 'lnLambda')

# ============================================================================
# Time stepping
# ============================================================================
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

print(f"\nTime Stepping:")
print(f"  t_max = {tMax} s")
print(f"  Nt = {Nt}")
print(f"  dt = {tMax/Nt:.2e} s")

# ============================================================================
# Output settings
# ============================================================================
ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('../outputs/quasilinear_whistler_output.h5')

# Save settings
settings_file = 'quasilinear_whistler_settings.h5'
ds.save(settings_file)

print("\n" + "="*70)
print("✓ Settings saved successfully!")
print("="*70)
print(f"\nSettings file: {settings_file}")
print(f"Output file: ../outputs/quasilinear_whistler_output.h5")
print("\nTo run manually:")
print(f"  /data/zhzhou/DREAM/build/iface/dreami {settings_file}")
print("="*70)

# Run the simulation
print("\nRunning simulation...")
do = runiface(ds, '../outputs/quasilinear_whistler_output.h5', quiet=False)

print("\n" + "="*70)
print("Quasilinear diffusion simulation completed!")
print("="*70)
