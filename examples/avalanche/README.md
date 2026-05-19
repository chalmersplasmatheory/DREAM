# Avalanche Runaway Electron Simulation

This directory contains a complete workflow for simulating avalanche runaway electron generation in DREAM with kinetic models.

## Directory Structure

```
avalanche/
├── scripts/           # Python scripts and shell script
│   ├── run_DIIID_fre.sh          # Main execution script
│   ├── generate_with_fre.py      # DREAM simulation setup and execution
│   ├── analyze_dreicer.py        # Basic analysis (density and growth rate)
│   ├── plot_dreicer_dist.py      # Plot distribution function
│   ├── plot_f_time_evolution.py  # Plot time evolution
│   ├── plot_p_perp_parallel.py   # Plot f_hot in p_perp-p_parallel space
│   ├── plot_fre_p_perp_parallel.py # Plot f_re in p_perp-p_parallel space
│   ├── visual.py                 # Comprehensive visualization
│   └── plot_logf_px_contour.py   # Plot log(f) contour in p-xi space
├── outputs/           # Simulation output files (HDF5)
└── figures/           # Generated plots and figures (PNG)
```

## Quick Start

1. **Run the complete workflow:**
   ```bash
   cd /data/zhzhou/DREAM/examples/avalanche/scripts
   ./run_DIIID_fre.sh
   ```

2. **Or run steps individually:**
   ```bash
   # Step 1: Generate settings and run simulation
   python generate_with_fre.py
   
   # Step 2: Analyze results
   python analyze_dreicer.py
   
   # Step 3: Visualize
   python plot_logf_px_contour.py
   # ... other plotting scripts
   ```

## Simulation Parameters

The simulation uses the following parameters (defined in `generate_with_fre.py`):

- **Electric field**: E = 0.05 V/m
- **Electron density**: n = 5×10¹⁸ m⁻³
- **Temperature**: T = 2165 eV (2.165 keV)
- **Magnetic field**: B = 1.4 T
- **Minor radius**: a = 0.67 m
- **Major radius**: R₀ = 1.67 m
- **Simulation time**: t_max = 100 ms
- **Time steps**: Nt = 5000

## Physics Models

- **Dreicer generation**: Neural network model (most accurate)
- **Avalanche generation**: Kinetic mode (most accurate, requires pCutAvalanche)
- **Hot electron grid**: Enabled (p_max = 1 m_e*c, Np = 100, Nxi = 40)
- **Runaway electron grid**: Enabled (p_max = 10 m_e*c, Np = 100, Nxi = 40)
- **Collision model**: Partially screened, full collision frequency

## Output Files

- `dreicer_with_fre_settings.h5` - DREAM input settings
- `dreicer_with_fre_output.h5` - Simulation results (in `outputs/`)
- Various PNG figures (in `figures/`)

## Visualization Scripts

| Script | Description |
|--------|-------------|
| `analyze_dreicer.py` | Runaway density and growth rate vs time |
| `plot_dreicer_dist.py` | Distribution function at specific time |
| `plot_f_time_evolution.py` | Time evolution of distribution |
| `plot_p_perp_parallel.py` | f_hot in perpendicular-parallel momentum space |
| `plot_fre_p_perp_parallel.py` | f_re in perpendicular-parallel momentum space |
| `plot_logf_px_contour.py` | **log(f) contour in p-xi space** (new!) |
| `visual.py` | Comprehensive multi-panel visualization |

## Notes

- The kinetic avalanche mode (`AVALANCHE_MODE_KINETIC`) is computationally intensive
- `pCutAvalanche = 2.0` sets the momentum threshold for avalanche participation
- Adjust grid resolution (Np, Nxi) and time steps (Nt) based on available computational resources
