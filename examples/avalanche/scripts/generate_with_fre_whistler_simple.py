#!/usr/bin/env python3
"""
Test script for quasilinear diffusion with simplified whistler dispersion relation.
Uses the fast analytical formula: ω = k|k_∥| * w
"""

import sys
import os
sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.CollisionHandler as Collisions
import numpy as np

# ============================================================================
# Physical parameters (matching QUADRE simulation)
# ============================================================================
E = 0.045 # Electric field strength (V/m), E/Ec = 0.045
n = 5e18           # Electron density (m^-3)
T = 2165           # Temperature (eV) - 2.165 keV
B0 = 1.4           # Magnetic field (Tesla)

# ============================================================================
# Create DREAM settings
# ============================================================================
ds = DREAMSettings()

# Set grid parameters
Nr = 2
Nxi = 20
Np = 40
pMax = 50

# Disable hot-tail grid (not needed for this test)
ds.hottailgrid.setEnabled(False)

# Set radial grid parameters
ds.radialgrid.setNr(2)  # Minimal radial resolution
ds.radialgrid.setMinorRadius(0.01)  # a = 1 cm
ds.radialgrid.setMajorRadius(1.67)  # R = 1.67 m
ds.radialgrid.setWallRadius(0.01)   # b = a (no scrape-off layer)
ds.radialgrid.setB0(B0)

# Enable runaway grid
ds.runawaygrid.setEnabled(True)
ds.runawaygrid.setNxi(Nxi)
ds.runawaygrid.setNp(Np)
ds.runawaygrid.setPmax(pMax)

# Set electric field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ion species (fully ionized deuterium)
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Set initial electron distribution
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

# Set collision mode
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_FULL
ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED 
ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS

# ============================================================================
# QUADRE wave parameters from simulation.log
# ============================================================================
quadre_wave_params = {
    'k_main': 54.58,              # Total wavenumber (m^-1)
    'ktheta_main': 2.42,          # Propagation angle (rad)
    'k_range': [50.61, 58.55],    # Wavenumber range (m^-1)
    'ktheta_range': [2.35, 2.49]  # Angle range (rad)
}

print(f"\nQUADRE Wave Parameters:")
print(f"  k_main = {quadre_wave_params['k_main']} m^-1")
print(f"  θ_k_main = {quadre_wave_params['ktheta_main']} rad ({quadre_wave_params['ktheta_main']*180/np.pi:.1f}°)")
print(f"  k range: [{quadre_wave_params['k_range'][0]:.2f}, {quadre_wave_params['k_range'][1]:.2f}] m^-1")
print(f"  θ_k range: [{quadre_wave_params['ktheta_range'][0]:.2f}, {quadre_wave_params['ktheta_range'][1]:.2f}] rad")

# ============================================================================
# Enable quasilinear diffusion with SIMPLIFIED dispersion relation
# ============================================================================
ds.eqsys.f_re.setQuasilinearDiffusion(
    enabled=True,
    quadre_params=quadre_wave_params,
    num_k=8,              # Match QUADRE discretization
    num_ktheta=20,        # Match QUADRE discretization
    amplitude=1e-10,      # Wave amplitude (normalized units)
    harmonic_mode='both',  # Include n = -2,-1,0,+1,+2
    use_simple_dispersion=True  # Use simplified dispersion: ω = k|k_∥| * w
)

print(f"  ✓ Quasilinear diffusion enabled (SIMPLIFIED dispersion)")
print(f"  - Spectrum: uniform grid")
print(f"  - Modes: {8} × {20} = {8*20}")
print(f"  - Harmonics: n ∈ {{-2, -1, 0, +1, +2}}")
print(f"  - Dispersion: ω = k|k_∥| * w (fast analytical formula)")

# ============================================================================
# Time stepping
# ============================================================================
t_max = 0.1  # seconds
Nt = 100
ds.timestep.setTmax(t_max)
ds.timestep.setNt(Nt)

# ============================================================================
# Save and run
# ============================================================================
output_file = '../outputs/quasilinear_whistler_simple.h5'
settings_file = 'quasilinear_whistler_simple_settings.h5'

ds.save(settings_file)

print(f"\n{'='*70}")
print(f"DREAM Quasilinear Diffusion Test - Simplified Dispersion Relation")
print(f"{'='*70}")
print(f"\nPlasma Parameters:")
print(f"  B0 = {B0} T")
print(f"  n_e = {n:.2e} m^-3")
print(f"  T = {T} eV")
print(f"  E = {E:.2f} V/m (E/Ec = 0.045)")
print(f"\nQuasilinear Diffusion Configuration:")
print(f"  Wave frequency: ~476 MHz (estimated from simplified dispersion)")
print(f"  k_parallel: -41 m^-1")
print(f"  k range: [{quadre_wave_params['k_range'][0]:.2f}, {quadre_wave_params['k_range'][1]:.2f}] m^-1")
print(f"  θ_k range: [{quadre_wave_params['ktheta_range'][0]:.2f}, {quadre_wave_params['ktheta_range'][1]:.2f}] rad")
print(f"  ✓ Quasilinear diffusion enabled (FAST MODE)")
print(f"  - Spectrum: uniform grid")
print(f"  - Modes: {8} × {20} = {8*20}")
print(f"  - Harmonics: n ∈ {{-2, -1, 0, +1, +2}}")
print(f"  - Dispersion: ω = k|k_∥| * w")
print(f"\nTime Stepping:")
print(f"  t_max = {t_max} s")
print(f"  Nt = {Nt}")
print(f"  dt = {t_max/Nt:.2e} s")
print(f"\n{'='*70}")
print(f"✓ Settings saved successfully!")
print(f"{'='*70}\n")
print(f"Settings file: {settings_file}")
print(f"Output file: {output_file}")
print(f"\nTo run manually:")
print(f"  /data/zhzhou/DREAM/build/iface/dreami {settings_file}")
print(f"{'='*70}\n")

# Run simulation
print("Running simulation...")
os.system(f'/data/zhzhou/DREAM/build/iface/dreami {settings_file}')
