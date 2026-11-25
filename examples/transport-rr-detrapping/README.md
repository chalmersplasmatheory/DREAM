# Isotropic RR detrapping validation

Run the example from this directory with

```bash
./run.py
```

The example compares ordinary Rechester–Rosenbluth (RR) diffusion with the
collisional trapped-pitch eigenmode in a shaped toroidal equilibrium. It runs
both models at `δB/B = 3e-4, 3e-3, 3e-2`. The default plasma has edge
elongation 1.55 and an initially neutral neon density of `1e19 m⁻³`. The neon
evolves with kinetic ionization and its approximate Jacobian.

The code reconstructs the trapped boundary from the analytic equilibrium and
computes its lowest pitch-scattering eigenvalue with finite differences.
Combining that eigenvalue with DREAM's radial-face scattering frequency gives
an independent prediction of the cap, `D_detrap = a² κ max(ν_D, 0) / j₀,₁²`.
The example checks that cap, the harmonic combination with ordinary RR, and
recovery of ordinary RR on the magnetic axis. Errors and the figure list are
written to `output/validation_summary.md`.

Figures in `output/figures/` show the coefficients at the outer face, the
transition between RR and detrapping limits, the predicted and measured cap,
and the volume-integrated hot-electron inventory. Use `--skip-plots` to run the
coefficient checks alone, or `--no-usetex` to render the figures without LaTeX.

For a time-dependent perturbation, run

```bash
./timescale_evolution.py
```

This example ramps `δB/B` through the transition and plots both the diffusion
coefficients and their cylindrical loss times `a²/(j₀,₁² D)`. The detrapping
time is `1/(κ ν_D)` on a trapped face. The script selects a radius and momentum
where the two limits cross when available. Use `--ir`, `--ip` or `--p` to select
the point manually. The plot is written to
`output_timescale/figures/eigenmode_timescale.pdf`.
