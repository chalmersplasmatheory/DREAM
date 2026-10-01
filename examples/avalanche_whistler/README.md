# Avalanche Runaway Electron Simulation with Whistler Wave Quasilinear Diffusion

This directory contains a complete workflow for simulating avalanche runaway electron generation in DREAM with kinetic models, **including quasilinear diffusion from externally injected whistler waves**.

## Directory Structure

```
avalanche_whistler/
├── scripts/                    # Python scripts and shell script
│   ├── run_DIIID_fre.sh                   # Main execution script (with/without waves)
│   ├── generate_with_fre.py               # Baseline simulation (no waves)
│   ├── generate_with_fre_whistler.py      # ⭐ Whistler wave simulation (NEW!)
│   ├── analyze_dreicer.py                 # Basic analysis (density and growth rate)
│   ├── plot_dreicer_dist.py               # Plot distribution function
│   ├── plot_f_time_evolution.py           # Plot time evolution
│   ├── plot_p_perp_parallel.py            # Plot f_hot in p_perp-p_parallel space
│   ├── plot_fre_p_perp_parallel.py        # Plot f_re in p_perp-p_parallel space
│   ├── visual.py                          # Comprehensive visualization
│   └── plot_logf_px_contour.py            # Plot log(f) contour in p-xi space
├── outputs/                    # Simulation output files (HDF5)
└── figures/                    # Generated plots and figures (PNG)
```

## Quick Start

### Option 1: Run with Whistler Waves (Recommended)

**Run the complete workflow with quasilinear diffusion:**
```bash
cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts
python generate_with_fre_whistler.py --amplitude 1e-5 --a 0.6 --R 1.67
```

This will:
1. Configure DREAM with 476 MHz whistler wave parameters (from QUADRE)
2. Enable quasilinear diffusion on runaway electrons
3. Run the simulation automatically
4. Save results to `../outputs/quasilinear_whistler_output.h5`

**Customize wave amplitude and geometry:**
```bash
# Test different wave strengths
python generate_with_fre_whistler.py --amplitude 1e-6    # Weaker waves
python generate_with_fre_whistler.py --amplitude 1e-4    # Stronger waves

# Adjust plasma geometry
python generate_with_fre_whistler.py --a 0.4 --R 1.67    # Smaller minor radius
```

### Option 2: Run Baseline (No Waves)

**For comparison, run without quasilinear diffusion:**
```bash
./run_DIIID_fre.sh
```

This runs the original `generate_with_fre.py` script (no waves) and generates all analysis plots.

### Option 3: Manual Step-by-Step Execution

```bash
# Step 1: Generate settings and run simulation (with waves)
python generate_with_fre_whistler.py --amplitude 1e-5

# Step 2: Analyze results
python analyze_dreicer.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures

# Step 3: Visualize
python plot_logf_px_contour.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures
# ... other plotting scripts
```

## Simulation Parameters

### Baseline Configuration (generate_with_fre.py)

The baseline simulation uses the following parameters:

- **Electric field**: E = 0.05 V/m
- **Electron density**: n = 5×10¹⁸ m⁻³
- **Temperature**: T = 2165 eV (2.165 keV)
- **Magnetic field**: B = 1.4 T
- **Minor radius**: a = 0.67 m
- **Major radius**: R₀ = 1.67 m
- **Simulation time**: t_max = 100 ms
- **Time steps**: Nt = 5000

### Whistler Wave Configuration (generate_with_fre_whistler.py)

**Wave Parameters (from QUADRE simulation):**
- **Wave frequency**: 476 MHz
- **Main wavenumber**: k = 54.58 m⁻¹
- **Propagation angle**: θ_k = 2.42 rad (obtuse angle for backward propagation)
- **k range**: [52.61, 56.55] m⁻¹
- **θ_k range**: [2.4, 2.45] rad
- **Wave spectrum**: Uniform grid, 8 × 20 = 160 modes
- **Harmonics**: n ∈ {-2, -1, 0, +1, +2}
- **Amplitude**: A = 1e-5 (normalized units, δB/B₀ ~ 10⁻⁵)

**Grid Parameters:**
- **Hot-tail grid**: p_max = 1 m_e*c, Np = 100, Nxi = 40
- **Runaway grid**: p_max = 50 m_e*c, Np = 200, Nxi = 40
- **Trapped-passing boundary layer**: dxi_max = 0.1
- **Simulation time**: t_max = 3.0 s
- **Time steps**: Nt = 3000
- **Radial points**: nr = 1 (single flux surface)

**Command-line Options:**
```bash
--amplitude AMPLITUDE   Wave amplitude (default: 1e-5)
--a MINOR_RADIUS         Minor radius in meters (default: 0.6)
--R MAJOR_RADIUS         Major radius in meters (default: 1.67)
```

## Physics Models

### Common Models (Both Configurations)

- **Dreicer generation**: Disabled (kinetic simulation captures it naturally via f_hot → f_re transfer)
- **Avalanche generation**: Kinetic mode (most accurate, requires pCutAvalanche = 2.0)
- **Collision model**: Partially screened, full collision frequency, energy-dependent Coulomb logarithm
- **Bremsstrahlung**: Stopping power mode
- **Synchrotron radiation**: Enabled for both f_hot and f_re

### Additional Physics (Whistler Wave Configuration)

- **Quasilinear diffusion**: ⭐ **Enabled on f_re**
  - Wave-particle resonance: ω - k_∥v_∥ + nΩ_ce/γ = 0
  - Resonant momenta: p_res ~ 17-60 m_e*c (depending on harmonic and pitch angle)
  - Diffusion tensor components: D_pp, D_pξ, D_ξξ (all non-zero)
  - Amplitude scaling: D(t) ∝ |A(t)|²
  - Pre-computed dispersion relation (simplified: ω = k|k_∥|·w)
  - Harmonic modes: n = -2, -1, 0, +1, +2

**Physical Mechanism:**
Whistler waves scatter runaway electrons in momentum space via cyclotron resonance, suppressing their growth by:
1. Pitch-angle scattering into loss cone
2. Momentum diffusion reducing acceleration efficiency
3. Energy redistribution across harmonics

## Output Files

### Settings Files
- `dreicer_with_fre_settings.h5` - Baseline input settings (no waves)
- `quasilinear_whistler_settings.h5` - Whistler wave input settings

### Simulation Results (in `outputs/`)
- `dreicer_with_fre_output.h5` - Baseline simulation results
- `quasilinear_whistler_output.h5` - Whistler wave simulation results

### Analysis Outputs (in `figures/`)
- Various PNG figures generated by plotting scripts
- Compare baseline vs. wave cases to see quasilinear diffusion effects

## Visualization Scripts

All scripts support both baseline and whistler wave simulations via `--data_dir` argument.

| Script | Description |
|--------|-------------|
| `analyze_dreicer.py` | Runaway density and growth rate vs time |
| `plot_dreicer_dist.py` | Distribution function at specific time |
| `plot_f_time_evolution.py` | Time evolution of distribution |
| `plot_p_perp_parallel.py` | f_hot in perpendicular-parallel momentum space |
| `plot_fre_p_perp_parallel.py` | f_re in perpendicular-parallel momentum space |
| `plot_logf_px_contour.py` | **log(f) contour in p-xi space** (recommended for QL diffusion!) |
| `visual.py` | Comprehensive multi-panel visualization |

**Example Usage:**
```bash
# Analyze whistler wave simulation
python analyze_dreicer.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures

# Plot log(f) contours to see diffusion effects
python plot_logf_px_contour.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures

# Compare with baseline
python plot_logf_px_contour.py --data_dir ../outputs/dreicer_with_fre_output.h5 --plot_dir ../figures
```

## Notes & Best Practices

### Numerical Stability

- **Wave amplitude tuning**: Start with small amplitudes (1e-6 to 1e-5) to avoid numerical instability
  - Too large amplitude → excessive diffusion → negative densities → crash
  - Too small amplitude → no visible effect on runaway electrons
  - Physical range: δB/B₀ ~ 10⁻⁵ (from Zehua Guo 2018 POP)

- **Time step selection**: dt = t_max/Nt should resolve the diffusion timescale
  - Check: dt < τ_diffusion ~ p²/D_pp
  - If simulation crashes with NaN, reduce dt or amplitude

### Grid Resolution

- **Momentum grid**: Np_re = 200 is recommended for resolving resonances at p ~ 17-60
- **Pitch-angle grid**: Nxi = 40 provides adequate resolution for trapped-passing dynamics
- **Trapped-passing boundary layer**: dxi_max = 0.1 ensures smooth transition

### Comparing With/Without Waves

To isolate the effect of quasilinear diffusion:

1. Run baseline: `python generate_with_fre.py`
2. Run with waves: `python generate_with_fre_whistler.py --amplitude 1e-5`
3. Compare outputs:
   ```bash
   python analyze_dreicer.py --data_dir ../outputs/dreicer_with_fre_output.h5
   python analyze_dreicer.py --data_dir ../outputs/quasilinear_whistler_output.h5
   ```

Look for:
- Reduced runaway electron density (n_re)
- Modified distribution function shape (especially near resonant momenta)
- Changes in growth rate over time

### Performance

- The kinetic avalanche mode (`AVALANCHE_MODE_KINETIC`) is computationally intensive
- Quasilinear diffusion adds minimal overhead when using pre-computed dispersion (~milliseconds per timestep)
- Total runtime: ~minutes for 3000 timesteps on a single CPU core

### Debugging

If you encounter issues:

1. **Check diffusion coefficients**:
   ```python
   import h5py
   with h5py.File('../outputs/quasilinear_whistler_output.h5', 'r') as f:
       # Look for validation output in terminal during simulation
       pass
   ```

2. **Verify resonance locations**:
   - Resonance should occur at p ~ 17-60 for the given wave parameters
   - Check debug output for "DEBUG findResonantP" messages

3. **Adjust parameters**:
   - Reduce amplitude if simulation crashes
   - Increase Np/Nxi if resonances are not well-resolved
   - Check that pMax_re > p_resonant for all harmonics
