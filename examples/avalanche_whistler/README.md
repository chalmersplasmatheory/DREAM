# Runaway Electron Simulation with Whistler Wave Quasilinear Diffusion

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
│   ├── plot_ql_diffusion_coefficients.py  # Quasilinear diffusion coefficient visualization
│   ├── dream_green_processor.py           # Green's function × f_re batch processor
│   ├── dream_se_power.py                  # Synchrotron radiation total power
|   ├── dump_hdf5_structure.py             # View quasilinear settings
|   ├── compare_dxixi_estimate.py          # Estimate diffusion time
├── outputs/                    # Simulation output files (HDF5)
└── figures/                    # Generated plots and figures (PNG)
```


### Run with Whistler Waves 

**Run the complete workflow with quasilinear diffusion:**
Change path to your directory! Bash command is also in run_DIIID_fre.sh 
```bash
cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts
PLOT_DIR="../figures_wk40_A4_t3_s2_i0.2_T1000"
DATA_DIR="${PLOT_DIR}/quasilinear_whistler_output.h5"
# Create plot directory if it doesn't exist
if [ ! -d "$PLOT_DIR" ]; then
    mkdir -p "$PLOT_DIR"
    echo "Created directory: $PLOT_DIR"
fi

python generate_with_fre_whistler.py \
    --amplitude 1e4 --a 0.02 --R 1.67 --T 100 \
    --tMax 2.5 --Nt 2500 --source off \
    --output $DATA_DIR \
    --start-inject-time 2.0 --inject-cycle-duration 0.2 --ramp-time 0.02 \
    --k-min 45 --k-max 60 \
    --ktheta-min 2.3 --ktheta-max 2.5 --Nxi 80 --num-k 40 --Np-hot 50 --Np-re 100
```
════════════════════════════════════════════════════════════════════════
Parameter quick reference
════════════════════════════════════════════════════════════════════════
(A) Quick command-line parameters — just edit the values in the call
    below; no need to touch the Python script:

  [Wave injection] (the most important ones)

    --amplitude              Wave amplitude (diffusion coefficient D ∝ A²)
    --start-inject-time      Wave injection start time (s)
    --inject-cycle-duration  ON+OFF injection cycle (s, 50% duty cycle)
    --ramp-time              Ramp-up time of each ON cycle (s)

  [Wave spectrum]

    --k-min / --k-max        Wavenumber range (m^-1)
    --ktheta-min / --ktheta-max   Propagation angle range (rad; obtuse = backward)
    --num-k / --num-ktheta   Spectral discretization (number of modes = product)

  [Plasma / geometry]

    --E electric field (V/m)   --n density (m^-3)   --T temperature (eV)
    --B0 magnetic field (T)    --a minor radius (m) --R major radius (m)

  [Grid / time]

    --Np-hot / --Np-re / --Nxi   Momentum and pitch grid points
    --tMax / --Nt                Simulation time (s) / number of time steps

  [Other]

    --source on|off          Kinetic avalanche source switch
    --output                 Output h5 path (recipe is saved next to it)

(B) Hard-coded values that require editing generate_with_fre_whistler.py
    (search for "HARDCODED" in the script to locate them):

    - pMax_hot=1 / pMax_re=50      Momentum grid limits (m_e c)
    - k_main=54.58 / ktheta_main=2.42   Central wavenumber / propagation angle
                                       (quadre_wave_params)
    - harmonic_mode='both'         Harmonic set n ∈ {-2..+2}
    - use_simple_dispersion=False  Full PDRK dispersion (True = simplified debug)
    - pCutAvalanche=2.0            Avalanche cutoff momentum
    - nr=1 and custom radial grid r_f   Number of radial points
    - setNumberOfSaveSteps(Nt//10) Number of saved output frames
    - Collision model / solver / boundary conditions / synchrotron switch, etc.

    NOTE: k_main/ktheta_main in quadre_wave_params are documentation only
    (never read by the settings API). The wave frequency (476 MHz) and
    k_par = -41 m^-1 are NOT inputs either: k_par = k*cos(theta_k) is derived
    per mode in C++, and the frequency follows from the whistler dispersion
    relation evaluated at (k, theta_k) together with B0 and n_e.
    The frequency is determined by the k spectrum.
════════════════════════════════════════════════════════════════════════

**Run the whistler-wave simulation followed by all post-processing plots:**
```bash
./run_DIIID_fre.sh
```

This runs `generate_with_fre_whistler.py`, then the consolidated plotting scripts
(`plot_f.py`, `plot_nre.py`, `plot_p_perp_parallel.py`). A commented block at the top
of the script shows how to run the baseline `generate_with_fre.py` instead.

### Visualization

```bash
python plot_f.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR \
    --species f_re \
    --t-min 1.0 --t-max 3.0 --dt 0.01 \
    --time_series 1.5,2.1 \
    --figures all

python plot_f.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR \
    --species f_hot --figures contour,spectra

# Runaway density, growth rate, normalized growth rate vs Rosenbluth 2*(E/Ec-1)
python plot_nre.py --data_dir $DATA_DIR --plot_dir $PLOT_DIR
```
