#!/usr/bin/env python3
"""Time-dependent diagnostic for the finite-detrapping RR limiter.

This script runs one eigenmode-limited Rechester-Rosenbluth case with a prescribed
time ramp in `dB/B`. It then reconstructs the two limiting coefficients

    D_RR(t)       = pi R0 (dB/B)^2 v G,
    D_detrap(t)   = a^2 kappa(xi_T) max(nu_D, 0) / j01^2,

where `kappa` is computed on an independent finite-difference pitch grid
and `nu_D` is the deflection frequency from DREAM output. The resulting figure shows both the coefficient
evolution and the equivalent loss timescales.
"""

import argparse
import os
from pathlib import Path
import sys
import warnings

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

MPLCONFIGDIR = Path(os.environ.get("TMPDIR", "/tmp")) / "dream-rr-detrapping-matplotlib"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(REPO_ROOT / "py"))

import numpy as np
from scipy.constants import c

from DREAM import DREAMOutput, DREAMSettings, runiface

import run as rr_example


def dBB_profile(t, args):
    """Logarithmic magnetic-perturbation ramp."""
    t = np.asarray(t, dtype=float)
    if args.dBB_min <= 0 or args.dBB_max <= 0:
        raise ValueError("dBB_min and dBB_max must both be positive.")

    if args.tmax <= 0:
        return args.dBB_min * np.ones_like(t)

    x = np.clip(t / args.tmax, 0.0, 1.0)
    return args.dBB_min * (args.dBB_max / args.dBB_min) ** x


def configure_time_dependent_transport(ds, args):
    r = np.linspace(0.0, args.minor_radius, args.nr + 1)
    t = np.linspace(0.0, args.tmax, args.ndbb)
    dBB = dBB_profile(t, args)[:, None] * np.ones((t.size, r.size))

    ds.eqsys.f_hot.transport.setDetrappingRechesterRosenbluth(
        dBB,
        t=t,
        r=r,
    )


def make_settings(args, output_path):
    ds = DREAMSettings()
    profiles = rr_example.initial_profiles(args)

    rr_example.configure_collisions(ds)
    rr_example.configure_radial_grid(ds, args)
    rr_example.configure_plasma(ds, profiles, args)
    rr_example.disable_runaway_sources(ds)
    rr_example.configure_hot_tail(ds, profiles, args)
    configure_time_dependent_transport(ds, args)
    rr_example.configure_solver_and_output(ds, output_path, args)

    ds.other.include("hottail/nu_D_f2")

    return ds


def run_case(args, outdir):
    settings_path = outdir / "eigenmode_timescale_settings.h5"
    output_path = outdir / "eigenmode_timescale_output.h5"

    ds = make_settings(args, output_path)
    ds.save(str(settings_path))

    if args.no_run:
        print(f"settings written: {settings_path}")
    else:
        print(f"running: {output_path}")
        runiface(ds, str(output_path), quiet=args.quiet)
        print(f"output written: {output_path}")

    return settings_path, output_path


def output_times(output, ntime):
    t_grid = np.asarray(output.grid.t, dtype=float)
    if t_grid.size == ntime:
        return t_grid
    if t_grid.size == ntime + 1:
        return t_grid[1:]
    return np.linspace(t_grid[0], t_grid[-1], ntime)


def extract_drr(output):
    Drr = np.asarray(output.other.hottail.Drr[:], dtype=float)
    if Drr.ndim == 4:
        if Drr.shape[2] != 1:
            raise RuntimeError(f"Expected isotropic Drr with one pitch cell, got {Drr.shape}.")
        return Drr[:, :, 0, :]
    if Drr.ndim == 3:
        return Drr
    raise RuntimeError(f"Unexpected Drr shape: {Drr.shape}.")


def extract_nuD_f2(output):
    if not hasattr(output.other.hottail, "nu_D_f2"):
        raise RuntimeError(
            "This output does not contain hottail/nu_D_f2. "
            "Run this script instead of reusing an older coefficient-only output."
        )

    nuD = np.asarray(output.other.hottail.nu_D_f2[:], dtype=float)
    if nuD.ndim == 4:
        return nuD[:, :, 0, :]
    if nuD.ndim == 3:
        return nuD
    raise RuntimeError(f"Unexpected nu_D_f2 shape: {nuD.shape}.")


def interpolate_nuD_to_radial_faces(nuD_cell, output):
    """Match the radial-face interpolation used in TrappingLimitedRRTransport."""
    r = np.asarray(output.grid.r, dtype=float)
    r_f = np.asarray(output.grid.r_f, dtype=float)
    nr = r.size

    nuD_face = np.empty((nuD_cell.shape[0], nr + 1, nuD_cell.shape[2]))
    nuD_face[:, 0, :] = nuD_cell[:, 0, :]
    nuD_face[:, nr, :] = nuD_cell[:, nr - 1, :]

    dr = np.diff(r_f)
    dr_f = np.diff(r)
    for ir in range(1, nr):
        w_right = 0.5 * dr[ir - 1] / dr_f[ir - 1]
        nuD_face[:, ir, :] = (
            w_right * nuD_cell[:, ir, :]
            + (1.0 - w_right) * nuD_cell[:, ir - 1, :]
        )

    return nuD_face


def reconstruct_coefficients(output, args):
    Deff = extract_drr(output)
    nuD_face = interpolate_nuD_to_radial_faces(extract_nuD_f2(output), output)
    t = output_times(output, Deff.shape[0])
    p = np.asarray(output.grid.hottail.p, dtype=float)

    G, xiT, Ddetrap_per_nuD = rr_example.transport_geometry_factors(np.asarray(output.grid.r_f), args)
    v = c * p / np.sqrt(1.0 + p**2)
    prescribed = output.settings.eqsys.f_hot.transport
    dBB_t = np.interp(t, prescribed.dBB_t, prescribed.dBB[:, 0])

    Drr = np.pi * args.major_radius * dBB_t[:, None, None] ** 2 * G[None, :, None] * v[None, None, :]
    Ddetrap = Ddetrap_per_nuD[None, :, None] * np.maximum(nuD_face, 0.0)
    Deff_analytic = rr_example.expected_eigenmode_rr(Drr, Ddetrap, xiT > 0)

    return {
        "t": t,
        "p": p,
        "D_eff": Deff,
        "D_eff_analytic": Deff_analytic,
        "D_RR": Drr,
        "D_detrap": Ddetrap,
        "G": G,
        "xiT": xiT,
    }


def choose_trace(data, args):
    nr_f = data["D_eff"].shape[1]
    np1 = data["D_eff"].shape[2]

    if args.ir is not None:
        ir = args.ir if args.ir >= 0 else nr_f + args.ir
    else:
        ir = None

    if args.ip is not None:
        ip = args.ip if args.ip >= 0 else np1 + args.ip
    elif args.p is not None:
        ip = int(np.argmin(np.abs(data["p"] - args.p)))
    else:
        ip = None

    if ir is not None and ip is not None:
        return ir, ip, "user-selected"

    ratio = np.divide(
        data["D_RR"],
        data["D_detrap"],
        out=np.full_like(data["D_RR"], np.nan),
        where=data["D_detrap"] > 0,
    )
    log_ratio = np.log10(ratio)
    closest = np.min(np.where(np.isfinite(log_ratio), np.abs(log_ratio), np.inf), axis=0)

    if ir is None and ip is None:
        crosses = np.any(log_ratio[:-1, :, :] * log_ratio[1:, :, :] <= 0, axis=0)
        if np.any(crosses):
            closest[~crosses] = np.inf
            ir, ip = np.unravel_index(np.argmin(closest), closest.shape)
            return int(ir), int(ip), "auto-crossing"
        ir, ip = np.unravel_index(np.argmin(closest), closest.shape)
    elif ir is None:
        ir = np.argmin(closest[:, ip])
    else:
        ip = np.argmin(closest[ir, :])

    return int(ir), int(ip), "closest-to-crossing"


def crossing_time(t, Drr, Ddetrap):
    ratio = np.divide(Drr, Ddetrap, out=np.full_like(Drr, np.nan), where=Ddetrap > 0)
    y = np.log(ratio)
    good = np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return None

    t = t[good]
    y = y[good]
    crossings = np.where(y[:-1] * y[1:] <= 0)[0]
    if crossings.size == 0:
        return None

    k = crossings[0]
    if y[k + 1] == y[k]:
        return t[k]
    return t[k] - y[k] * (t[k + 1] - t[k]) / (y[k + 1] - y[k])


def dominant_label(Drr, Ddetrap):
    if Drr < Ddetrap:
        return "RR-limited"
    return "detrapping-limited"


def plot_timescale_evolution(data, ir, ip, source, args, figdir):
    import matplotlib.pyplot as plt

    t_us = 1e6 * data["t"]
    p = data["p"][ip]

    Drr = data["D_RR"][:, ir, ip]
    Ddetrap = data["D_detrap"][:, ir, ip]
    Deff = data["D_eff"][:, ir, ip]
    Deff_analytic = data["D_eff_analytic"][:, ir, ip]

    length_squared = args.minor_radius**2 / rr_example.J01_SQUARED
    tau_rr = length_squared / Drr
    tau_detrap = length_squared / Ddetrap
    tau_eff = length_squared / Deff

    tcross = crossing_time(data["t"], Drr, Ddetrap)

    fig, axs = plt.subplots(2, 1, figsize=(5.3, 5.0), sharex=True)

    ax = axs[0]
    ax.semilogy(t_us, Drr, color="#5b2a86", label=r"$D_{\rm RR}$ analytic")
    ax.semilogy(t_us, Ddetrap, color="#8c1d18", label=r"$D_{\rm detrap}$ analytic")
    ax.semilogy(t_us, Deff, color="#0f766e", label=r"$D_{\rm eff}$ DREAM")
    ax.semilogy(t_us, Deff_analytic, "--", color="#0f766e", label=r"$D_{\rm eff}$ harmonic")
    ax.set_ylabel(r"$D\;[\mathrm{m^2\,s^{-1}}]$")
    ax.legend(frameon=False, fontsize=7)

    ax = axs[1]
    ax.semilogy(t_us, tau_rr, color="#5b2a86", label=r"$\tau_{\rm RR}=a^2/(j_{0,1}^2D_{\rm RR})$")
    ax.semilogy(t_us, tau_detrap, color="#8c1d18", label=r"$\tau_{\rm detrap}=1/(\kappa\nu_D)$")
    ax.semilogy(t_us, tau_eff, color="#0f766e", label=r"$\tau_{\rm eff}=a^2/(j_{0,1}^2D_{\rm eff})$")
    ax.set_xlabel(r"$t\;[\mu\mathrm{s}]$")
    ax.set_ylabel(r"$\tau\;[\mathrm{s}]$")
    ax.legend(frameon=False, fontsize=7)

    if tcross is not None:
        for ax in axs:
            ax.axvline(1e6 * tcross, color="0.3", ls=":", lw=1.0)

    title = rf"$r$-face {ir}, $p={p:.3g}$, {source}"
    axs[0].set_title(title)

    path = figdir / "eigenmode_timescale.pdf"
    fig.savefig(path)
    plt.close(fig)

    return path, tcross


def print_summary(data, ir, ip, source, tcross):
    Drr = data["D_RR"][:, ir, ip]
    Ddetrap = data["D_detrap"][:, ir, ip]
    ratio = Drr / Ddetrap

    print("\nTimescale diagnostic")
    print(f"  trace selection     : {source}")
    print(f"  radial face index   : {ir}")
    print(f"  momentum index      : {ip}")
    print(f"  p/(m_e c)           : {data['p'][ip]:.6g}")
    print(f"  initial D_RR/D_detrap: {ratio[0]:.6e} ({dominant_label(Drr[0], Ddetrap[0])})")
    print(f"  final D_RR/D_detrap  : {ratio[-1]:.6e} ({dominant_label(Drr[-1], Ddetrap[-1])})")
    if tcross is None:
        print("  inversion time      : no crossing found")
    else:
        print(f"  inversion time      : {1e6 * tcross:.6g} us")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run and plot a time-dependent finite-detrapping RR timescale diagnostic.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--output-dir", default=str(SCRIPT_DIR / "output_timescale"),
        help="Directory for settings, output and figure.")
    parser.add_argument("--dBB-min", type=float, default=3.0e-4,
        help="Initial magnetic perturbation level.")
    parser.add_argument("--dBB-max", type=float, default=3.0e-2,
        help="Final magnetic perturbation level.")
    parser.add_argument("--ndbb", type=int, default=128,
        help="Number of samples used to prescribe the dB/B ramp.")

    parser.add_argument("--nr", type=int, default=10, help="Number of radial cells.")
    parser.add_argument("--np", type=int, default=90, help="Number of momentum cells.")
    parser.add_argument("--pmax", type=float, default=2.5, help="Maximum hot-tail momentum.")
    parser.add_argument("--tmax", type=float, default=3.0e-4, help="Simulation time [s].")
    parser.add_argument("--nt", type=int, default=80, help="Number of time steps.")

    parser.add_argument("--minor-radius", type=float, default=0.5, help="Minor radius [m].")
    parser.add_argument("--major-radius", type=float, default=3.0, help="Major radius [m].")
    parser.add_argument("--B0", type=float, default=5.0, help="On-axis toroidal magnetic field [T].")
    parser.add_argument("--plasma-current", type=float, default=8.0e5,
        help="Reference plasma current for the analytic poloidal flux [A].")
    parser.add_argument("--elongation", type=float, default=1.55,
        help="Edge elongation of the analytic toroidal flux surfaces.")
    parser.add_argument("--triangularity", type=float, default=0.0,
        help="Edge triangularity of the analytic toroidal flux surfaces.")
    parser.add_argument("--shafranov-shift", type=float, default=0.0,
        help="Edge Shafranov shift of the analytic toroidal flux surfaces [m].")

    parser.add_argument("--density", type=float, default=5.0e19, help="Cold electron density [m^-3].")
    parser.add_argument("--neon-density", type=float, default=1.0e19,
        help="Initial flat neutral neon density [m^-3]. Set to 0 to disable neon.")
    parser.add_argument("--cold-temperature", type=float, default=1.0e3, help="Cold electron temperature [eV].")
    parser.add_argument("--hot-temperature", type=float, default=7.5e4, help="Initial f_hot temperature [eV].")
    parser.add_argument("--hot-fraction", type=float, default=5.0e-3, help="Initial hot density fraction.")
    parser.add_argument("--hollowing", type=float, default=0.65, help="Quadratic radial hollowing of n_hot.")
    parser.add_argument("--electric-field", type=float, default=0.0, help="Prescribed electric field [V/m].")

    parser.add_argument("--ir", type=int, default=None,
        help="Radial face index to plot. Negative indices are accepted.")
    parser.add_argument("--ip", type=int, default=None,
        help="Momentum index to plot. Negative indices are accepted.")
    parser.add_argument("--p", type=float, default=None,
        help="Momentum value to plot if --ip is not given.")

    parser.add_argument("--no-run", action="store_true", help="Only write settings.")
    parser.add_argument("--plot-only", action="store_true", help="Use an existing output file.")
    parser.add_argument("--quiet", action="store_true", help="Run DREAM quietly.")
    parser.add_argument("--no-usetex", action="store_false", dest="usetex",
        help="Disable external LaTeX text rendering.")
    parser.set_defaults(usetex=True)

    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir).resolve()
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    output_path = outdir / "eigenmode_timescale_output.h5"
    if not args.plot_only:
        _, output_path = run_case(args, outdir)

    if args.no_run:
        return 0
    if not output_path.exists():
        raise SystemExit(f"Missing output file: {output_path}")

    rr_example.configure_matplotlib(usetex=args.usetex)
    output = DREAMOutput(str(output_path))
    data = reconstruct_coefficients(output, args)
    ir, ip, source = choose_trace(data, args)
    figure_path, tcross = plot_timescale_evolution(data, ir, ip, source, args, figdir)

    print_summary(data, ir, ip, source, tcross)
    print(f"  figure              : {figure_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
