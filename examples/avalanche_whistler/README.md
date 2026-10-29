# Avalanche Runaway Electron Simulation with Whistler Wave Quasilinear Diffusion

This directory contains a complete workflow for simulating avalanche runaway electron generation in DREAM with kinetic models, **including quasilinear diffusion from externally injected whistler waves**.

## Directory Structure

```
avalanche_whistler/
├── scripts/                    # Python scripts and shell scripts
│   ├── run_DIIID_fre.sh                   # Main execution script (simulation + post-processing)
│   ├── run_dream_green.sh                 # Synchrotron radiation pipeline (SOFT Green's function)
│   ├── generate_with_fre.py               # Baseline simulation (no waves)
│   ├── generate_with_fre_whistler.py      # ⭐ Whistler wave simulation
│   ├── plot_f.py                          # Distribution functions (p-xi contours + 1D spectra)
│   ├── plot_nre.py                        # Runaway density & growth rates
│   ├── plot_p_perp_parallel.py            # f_hot / f_re in p_perp-p_parallel space
│   ├── plot_ql_diffusion_coefficients.py  # Quasilinear diffusion coefficient visualization
│   ├── dream_green_processor.py           # Green's function × f_re batch processor
│   ├── dream_se_power.py                  # Synchrotron radiation total power
│   └── deprecated/                        # Legacy plotting scripts (kept for reference)
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

### Option 2: Run Full Workflow (SLURM)

**Run the whistler-wave simulation followed by all post-processing plots:**
```bash
./run_DIIID_fre.sh
```

This runs `generate_with_fre_whistler.py`, then the consolidated plotting scripts
(`plot_f.py`, `plot_nre.py`, `plot_p_perp_parallel.py`). A commented block at the top
of the script shows how to run the baseline `generate_with_fre.py` instead.

### Option 3: Manual Step-by-Step Execution

```bash
# Step 1: Generate settings and run simulation (with waves)
python generate_with_fre_whistler.py --amplitude 1e-5

# Step 2: Runaway density & growth rate analysis
python plot_nre.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures

# Step 3: Distribution functions (p-xi contours + 1D spectra)
python plot_f.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures --figures all

# Step 4: f_hot / f_re in (p_parallel, p_perp) space
python plot_p_perp_parallel.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures
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

**Key wave-injection parameters (as set in run_DIIID_fre.sh), by importance:**

1. **Wave strength** — the single most important knob (diffusion coefficient D(t) ∝ |A(t)|²):
   ```bash
   --amplitude 1e7
   ```
2. **Injection timing program** (must lie within `--tMax`):
   ```bash
   --start-inject-time 1.6        # when wave injection starts (s)
   --inject-cycle-duration 0.2    # ON+OFF cycle period (s)
   --ramp-time 0.02               # ramp-up smoothing to avoid step shocks (s)
   ```
3. **Wave spectrum**:
   ```bash
   --k-min 40 --k-max 60                # wavenumber range (m^-1)
   --ktheta-min 2.3 --ktheta-max 2.5    # propagation angle range (rad, obtuse → backward)
   --num-k 40 --num-ktheta 20           # spectral discretization (modes = product)
   ```
4. **Momentum-space resolution** (pitch grid resolves the resonance curves):
   ```bash
   --Nxi 80 --Np-re 100 --Np-hot 50
   ```
5. **Plasma parameters feeding the wave module**:
   ```bash
   --n 5e18    # density, passed as `density=` to the quasilinear module
   --B0 1.4    # magnetic field → cyclotron resonance location
   --E 0.05    # electric field (runaway drive being suppressed)
   ```

**Hard-coded wave parameters** (not exposed as command-line options):
- Wave frequency: 476 MHz
- Parallel wavenumber: k∥ = −41 m⁻¹
- Main wavenumber / angle: k_main = 54.58 m⁻¹, θ_main = 2.42 rad
- Harmonics: n ∈ {−2, −1, 0, +1, +2} (`harmonic_mode='both'`)
- Dispersion relation: full (`use_simple_dispersion=False`)
- Diffusion target: f_re (runaway distribution)
- Momentum grid limits: pMax_hot = 1, pMax_re = 50 (m_e c)

> **Note:** wave-injection parameters are **not stored in the output HDF5 file**
> (`settings/eqsys/f_re` contains no quasilinear entries and `eqsys/other` has no wave
> diagnostics — see `hdf5_structure.txt`). They are only recorded in the generator's
> command-line invocation and in `quasilinear_whistler_settings.h5`. Keep this in mind
> when analyzing old outputs.

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

Each generator also saves a copy of its settings file (the "recipe") **next to the
output HDF5 file**, so every run directory contains both its recipe and its results.

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
| `plot_f.py` | Distribution functions: log(f) contours in p-ξ space (with whistler resonance lines), 1D spectra f(p), pitch-angle comparison, trapped-electron spectrum |
| `plot_nre.py` | Runaway density (incl. trapped/passing decomposition), growth rate γ(t), normalized growth rate vs Rosenbluth 2(E/E_c−1); multiple input files give an E/E_c scan |
| `plot_p_perp_parallel.py` | f_hot and/or f_re in perpendicular-parallel momentum space (scatter, contour, 3-time evolution) |
| `plot_ql_diffusion_coefficients.py` | Quasilinear diffusion coefficient visualization |
| `dream_green_processor.py` + `dream_se_power.py` | Synchrotron radiation pipeline (Green's function × f_re → total power) |

`plot_f.py` supports selecting the data time range and the plotting time interval
(e.g. one figure every 0.01 s):

```bash
# All distribution figures at default times
python plot_f.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures --figures all

# One p-xi contour every 0.01 s between t = 1.0 and 2.0 s
python plot_f.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures \
    --t-min 1.0 --t-max 2.0 --dt 0.01 --figures contour

# Explicit time list for multi-time figures (seconds)
python plot_f.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures \
    --time_series 1.0,1.5,2.0 --figures spectra,grid
```

**Example Usage:**
```bash
# Runaway density and growth rate
python plot_nre.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures

# Distribution function contours to see wave-particle diffusion effects
python plot_f.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures --figures all

# f_hot and f_re in (p_parallel, p_perp) space
python plot_p_perp_parallel.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures --species both

# Multi-run E/Ec scan (several output files)
python plot_nre.py --data_dir out_E1.0.h5 out_E1.5.h5 out_E2.0.h5 --plot_dir ../figures_scan

# Compare with baseline
python plot_f.py --data_dir ../outputs/dreicer_with_fre_output.h5 --plot_dir ../figures
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
   python plot_nre.py --data_dir ../outputs/dreicer_with_fre_output.h5 --plot_dir ../figures
   python plot_nre.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures
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
