#!/usr/bin/env python3
#
# Runaway electron density and growth rate analysis (consolidated).
#
# Figures:
#   density_evolution.png       - n_cold/n_hot/n_re + trapped/passing + n_total (single-file)
#   nre_growth_rate.png         - n_re(t) and gamma(t) panels (single-file)
#   normalized_growth_rate.png  - gamma_norm(t) vs Rosenbluth 2*(E/Ec-1)
#   nre_vs_time_scan.png        - n_re(t) overlaid for several runs (multi-file)
#   growth_rate_vs_time_scan.png - gamma(t) overlaid (multi-file)
#   avg_growth_vs_Ec.png        - average gamma vs E/Ec scan (multi-file)
#
# Single-file mode:  one --data_dir h5 → detailed density/growth figures.
# Multi-run mode:    several --data_dir h5 files → E/Ec scan overlay figures.
#
# Merges visual.py + analyze_dreicer.py + plot_dreicer_dist.py (figs 2-3).
#
# Examples:
#   python plot_nre.py --data_dir ../outputs/quasilinear_whistler_output.h5 --plot_dir ../figures_x
#   python plot_nre.py --data_dir out_E1.0.h5 out_E1.5.h5 out_E2.0.h5 --plot_dir ../figures_scan
#

import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import DREAMOutput


def rosenbluth_norm(ef):
    """Rosenbluth-Parail normalized growth rate: gamma_norm = 2*(E/Ec - 1)."""
    return 2.0 * (ef - 1.0)


def extract_E_over_Ec(do):
    """Return (E/Ec value or None, legend label) for one DREAM output."""
    try:
        E0 = float(do.eqsys.E_field.get()[0, 0])
        Ec0 = float(do.other.fluid.Ectot[:][0, 0])
        ef = E0 / Ec0
        return ef, f'E/Ec = {ef:.2f}'
    except Exception as e:
        print(f"Warning: could not extract E/Ec ({e}); using filename as label.")
        return None, None


def load_runs(data_paths):
    """Load DREAMOutput objects for the given paths.

    Returns list of dicts: {'path', 'do', 'ef', 'label'} in CLI order.
    """
    runs = []
    for path in data_paths:
        if not os.path.exists(path):
            print(f"Warning: file not found, skipping: {path}")
            continue
        do = DREAMOutput(path)
        ef, ef_label = extract_E_over_Ec(do)
        label = ef_label if ef_label else Path(path).stem
        runs.append({'path': path, 'do': do, 'ef': ef, 'label': label})
        print(f"Loaded {path} → {label}")
    return runs


def _align_mid(x, nt, t_grid):
    """Align a per-time quantity to the nt-1 midpoint samples used for gamma."""
    x = np.asarray(x)
    nmid = nt - 1
    if x.shape[0] == nmid:
        return x
    if x.shape[0] == nt:
        return 0.5 * (x[:-1] + x[1:])
    t_x = np.linspace(t_grid[0], t_grid[-1], x.shape[0])
    t_avg = 0.5 * (t_grid[:-1] + t_grid[1:])
    return np.interp(t_avg, t_x, x)


def compute_growth_rates(do, n_re_floor=1e5):
    """
    Growth rate from n_re(t): gamma = (1/n_re) dn_re/dt, and the normalized
    form gamma_norm = gamma * tau_EERel * 2 * lnLambdaC (mentor convention,
    comparable with Rosenbluth 2*(E/Ec-1)).

    Returns (t_avg, gamma, gamma_norm) with invalid samples removed.
    """
    t = do.grid.t[:]
    n_re = do.eqsys.n_re.get()[:, 0]

    dt = np.diff(t)
    dn_re = np.diff(n_re)
    n_re_avg = 0.5 * (n_re[:-1] + n_re[1:])
    t_avg = 0.5 * (t[:-1] + t[1:])

    with np.errstate(divide='ignore', invalid='ignore'):
        gamma = dn_re / dt / n_re_avg

    mask = (n_re_avg > n_re_floor) & np.isfinite(gamma)
    if not np.any(mask):
        return np.array([]), np.array([]), np.array([])

    # Normalization factors (stored on their own time grid → align to midpoints)
    try:
        tau = _align_mid(do.other.fluid.tauEERel[:][:, 0], len(t), t)
        lnL = _align_mid(do.other.fluid.lnLambdaC[:][:, 0], len(t), t)
        gamma_norm = gamma * tau * 2.0 * lnL
    except Exception as e:
        print(f"Warning: could not compute normalized growth rate ({e})")
        gamma_norm = np.full_like(gamma, np.nan)

    return t_avg[mask], gamma[mask], gamma_norm[mask]


def compute_trapped_passing(do):
    """
    Trapped/passing decomposition of the runaway density from f_re
    (V'-weighted integration, rescaled so trapped+passing == fluid n_re).

    Returns (t, n_re, n_re_trapped, n_re_passing); trapped/passing are NaN
    when f_re is not available.
    """
    t = do.grid.t[:]
    n_re = do.eqsys.n_re.get()[:, 0]

    if not (hasattr(do.eqsys, 'f_re') and do.eqsys.f_re is not None):
        print("Warning: f_re not available; skipping trapped/passing decomposition.")
        return t, n_re, np.full_like(n_re, np.nan), np.full_like(n_re, np.nan)

    f_re_r0 = do.eqsys.f_re[:, 0, :, :]        # [nt, nxi, np]
    p = do.grid.runaway.p[:]
    xi = do.grid.runaway.xi[:]
    Vp = do.grid.runaway.Vprime[:]             # [nr, nxi, np]

    if not (hasattr(do.grid, 'xi0TrappedBoundary') and do.grid.xi0TrappedBoundary is not None):
        print("Warning: xi0TrappedBoundary not available; skipping trapped/passing split.")
        return t, n_re, np.full_like(n_re, np.nan), np.full_like(n_re, np.nan)

    xi0 = do.grid.xi0TrappedBoundary[0]
    mask_trapped = np.abs(xi) < xi0
    mask_passing = ~mask_trapped

    Vp_r0 = Vp[0, :, :]                        # [nxi, np]
    integrand = f_re_r0 * Vp_r0[np.newaxis, :, :]

    n_re_trapped = np.trapz(
        np.trapz(integrand[:, mask_trapped, :], p, axis=-1),
        xi[mask_trapped], axis=-1)
    n_re_passing = np.trapz(
        np.trapz(integrand[:, mask_passing, :], p, axis=-1),
        xi[mask_passing], axis=-1)

    total_from_f = n_re_trapped + n_re_passing
    scale = np.where(total_from_f > 0, n_re / total_from_f, 1.0)
    n_re_trapped = n_re_trapped * scale
    n_re_passing = n_re_passing * scale

    return t, n_re, n_re_trapped, n_re_passing


# ───────────────────────── figures (single-file) ─────────────────────────

def plot_density_evolution(do, plot_dir, t_min=0.5):
    """
    Density evolution: n_cold, n_hot, n_re, trapped/passing decomposition
    and n_total. Requires f_re for the trapped/passing split (skipped otherwise).
    """
    t, n_re, n_re_trapped, n_re_passing = compute_trapped_passing(do)
    n_cold = do.eqsys.n_cold[:, 0]
    n_hot = do.eqsys.n_hot[:, 0]

    mask = t >= t_min
    if not np.any(mask):
        mask = np.ones_like(t, dtype=bool)

    plt.figure(figsize=(7, 6))
    plt.plot(t[mask], n_cold[mask], label=r'$n_{\rm cold}$')
    plt.plot(t[mask], n_hot[mask], label=r'$n_{\rm hot}$')
    plt.plot(t[mask], n_re[mask], label=r'$n_{\rm re}$')
    if np.any(np.isfinite(n_re_trapped)):
        plt.plot(t[mask], n_re_trapped[mask], '--', label=r'$n_{\rm re, trapped}$')
        plt.plot(t[mask], n_re_passing[mask], ':', label=r'$n_{\rm re, passing}$')
    plt.plot(t[mask], (n_cold + n_hot + n_re)[mask], '-.', label=r'$n_{\rm total}$')

    plt.xlabel('Time (s)')
    plt.ylabel(r'Electron density (m$^{-3}$)')
    plt.title(f'Electron density evolution (r = 0, $t \\geq {t_min}$ s)\n'
              r'$n_{\rm cold}$, $n_{\rm hot}$, $n_{\rm re}$ and trapped/passing decomposition')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.yscale('log')
    plt.tight_layout()

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, 'density_evolution.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close()
    return save_path


def plot_nre_growth_2panel(do, plot_dir, t_min=0.5, n_re_floor=1e5):
    """Two panels: n_re(t) (log) and growth rate gamma(t) with final-value annotation."""
    t = do.grid.t[:]
    n_re = do.eqsys.n_re.get()[:, 0]

    t_avg, gamma, _ = compute_growth_rates(do, n_re_floor)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    mask = t >= t_min
    if not np.any(mask):
        mask = np.ones_like(t, dtype=bool)

    ax1.plot(t[mask], n_re[mask], 'b-', linewidth=2)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel(r'Runaway electron density $n_{\rm re}$ (m$^{-3}$)')
    ax1.set_title(f'Runaway electron density $n_{{\\rm re}}(t)$ ($t \\geq {t_min}$ s)')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_yscale('log')

    if len(t_avg) > 0:
        gmask = t_avg >= t_min
        if not np.any(gmask):
            gmask = np.ones_like(t_avg, dtype=bool)
        ax2.plot(t_avg[gmask], gamma[gmask], 'r-', linewidth=2)

        final_gamma = gamma[gmask][-1]
        final_time = t_avg[gmask][-1]
        ax2.annotate(f'Final $\\gamma$ = {final_gamma:.3e} s$^{{-1}}$\n'
                     f'at t = {final_time:.3f} s',
                     xy=(final_time, final_gamma),
                     xytext=(0.6 * final_time, 0.7 * final_gamma),
                     fontsize=11,
                     bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                     arrowprops=dict(arrowstyle='->',
                                     connectionstyle='arc3,rad=0.3', lw=2))
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel(r'Growth rate $\gamma = (1/n_{\rm re})\,dn_{\rm re}/dt$ (s$^{-1}$)')
    ax2.set_title(f'Runaway growth rate $\\gamma(t)$ ($t \\geq {t_min}$ s)')
    ax2.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, 'nre_growth_rate.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close(fig)
    return save_path


def plot_normalized_growth_rate(runs, plot_dir, t_min=0.5, n_re_floor=1e5):
    """gamma_norm(t) for 1..N runs with the Rosenbluth 2*(E/Ec-1) reference lines."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    plotted = False
    for run in runs:
        t_avg, _, gamma_norm = compute_growth_rates(run['do'], n_re_floor)
        if len(t_avg) == 0:
            print(f"Warning: no valid growth rate data for {run['label']}")
            continue
        mask_t = t_avg >= t_min
        if not np.any(mask_t):
            mask_t = np.ones_like(t_avg, dtype=bool)
        ax.plot(t_avg[mask_t], gamma_norm[mask_t], '-', linewidth=2,
                label=f"DREAM ({run['label']})")
        plotted = True

        if run['ef'] is not None:
            rp = rosenbluth_norm(run['ef'])
            ax.axhline(y=rp, linestyle='--', linewidth=1.5,
                       label=f"Rosenbluth $2(E/E_c-1)$ = {rp:.2f}")
            print(f"  {run['label']}: Rosenbluth 2*(E/Ec-1) = 2*({run['ef']:.4f}-1) = {rp:.4f}")

    if not plotted:
        plt.close(fig)
        return None

    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r'Normalized growth rate $\gamma\,\tau_{EE}\cdot 2\ln\Lambda_c$')
    ax.set_title('Normalized runaway growth rate vs time\n'
                 r'compared with Rosenbluth theory $\gamma_{\rm norm} = 2(E/E_c - 1)$')
    ax.legend(fontsize=11)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_xlim(left=t_min)
    plt.tight_layout()

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, 'normalized_growth_rate.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close(fig)
    return save_path


# ───────────────────────── figures (multi-run scan) ─────────────────────────

def plot_scan_overlays(runs, plot_dir, t_min=0.5, n_re_floor=1e5):
    """n_re(t) and gamma(t) overlaid for several runs (legend = E/Ec)."""
    saved = []

    # n_re(t)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    for run in runs:
        t = run['do'].grid.t[:]
        n_re = run['do'].eqsys.n_re.get()[:, 0]
        mask = t >= t_min
        if not np.any(mask):
            mask = np.ones_like(t, dtype=bool)
        ax.plot(t[mask], n_re[mask], linewidth=2, label=run['label'])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r'Runaway electron density $n_{\rm re}$ (m$^{-3}$)')
    ax.set_title(f'Runaway electron density $n_{{\\rm re}}(t)$ vs time for different $E/E_c$ '
                 f'($t \\geq {t_min}$ s)')
    ax.set_yscale('log')
    ax.legend(fontsize=11)
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    os.makedirs(plot_dir, exist_ok=True)
    p1 = os.path.join(plot_dir, 'nre_vs_time_scan.png')
    plt.savefig(p1, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {p1}")
    plt.close(fig)
    saved.append(p1)

    # gamma(t)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    for run in runs:
        t_avg, gamma, _ = compute_growth_rates(run['do'], n_re_floor)
        if len(t_avg) == 0:
            continue
        mask = t_avg >= t_min
        if not np.any(mask):
            mask = np.ones_like(t_avg, dtype=bool)
        ax.plot(t_avg[mask], gamma[mask], linewidth=2, label=run['label'])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r'Growth rate $\gamma = (1/n_{\rm re})\,dn_{\rm re}/dt$ (s$^{-1}$)')
    ax.set_title(f'Runaway growth rate $\\gamma(t)$ vs time for different $E/E_c$ '
                 f'($t \\geq {t_min}$ s)')
    ax.legend(fontsize=11)
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    p2 = os.path.join(plot_dir, 'growth_rate_vs_time_scan.png')
    plt.savefig(p2, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {p2}")
    plt.close(fig)
    saved.append(p2)

    return saved


def plot_avg_growth_vs_Ec(runs, plot_dir, n_re_floor=1e5):
    """Average growth rate vs E/Ec: physical gamma and normalized vs Rosenbluth."""
    ef_list, avg_gamma_list, avg_norm_list = [], [], []

    for run in runs:
        if run['ef'] is None:
            print(f"Warning: {run['label']} has no E/Ec value; skipping in scan plot.")
            continue
        t_avg, gamma, gamma_norm = compute_growth_rates(run['do'], n_re_floor)
        if len(gamma) < 10:
            print(f"Warning: too few valid growth-rate samples for {run['label']}; skipping.")
            continue
        avg_gamma = np.mean(gamma[-20:])
        avg_norm = np.nanmean(gamma_norm[-20:])
        ef_list.append(run['ef'])
        avg_gamma_list.append(avg_gamma)
        avg_norm_list.append(avg_norm)
        print(f"  {run['label']}: <gamma> = {avg_gamma:.4e} s^-1, "
              f"<gamma_norm> = {avg_norm:.4f}")

    if len(ef_list) == 0:
        print("Error: no runs with usable E/Ec and growth data; skipping scan figure.")
        return None

    order = np.argsort(ef_list)
    ef_arr = np.asarray(ef_list)[order]
    avg_gamma_arr = np.asarray(avg_gamma_list)[order]
    avg_norm_arr = np.asarray(avg_norm_list)[order]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.plot(ef_arr, avg_gamma_arr, 'bo-', markersize=8, label='DREAM simulation')
    ax1.set_xlabel('$E/E_c$')
    ax1.set_ylabel(r'Average growth rate $\langle\gamma\rangle$ (s$^{-1}$)')
    ax1.set_title('Average runaway growth rate vs electric field')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)

    ax2.plot(ef_arr, avg_norm_arr, 'bo-', markersize=8,
             label='DREAM (normalized)')
    ef_range = np.linspace(min(ef_arr), max(ef_arr), 100)
    ax2.plot(ef_range, rosenbluth_norm(ef_range), 'r--', linewidth=2,
             label=r'Rosenbluth: $2(E/E_c - 1)$')
    for ef, nr in zip(ef_arr, avg_norm_arr):
        ax2.annotate(f'{nr:.2f}', (ef, nr), textcoords='offset points',
                     xytext=(0, 10), ha='center', fontsize=9)
    ax2.set_xlabel('$E/E_c$')
    ax2.set_ylabel(r'Normalized growth rate $\gamma\,\tau_{EE}\cdot 2\ln\Lambda_c$')
    ax2.set_title('Normalized growth rate vs Rosenbluth theory\n'
                  r'$\gamma_{\rm norm} = 2(E/E_c - 1)$')
    ax2.legend(fontsize=11)
    ax2.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, 'avg_growth_vs_Ec.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close(fig)
    return save_path


# ───────────────────────── statistics ─────────────────────────

def print_stats(runs, n_re_floor=1e5, t_min=0.5):
    """Print growth-rate statistics for each run."""
    print("\n" + "=" * 70)
    print("GROWTH RATE STATISTICS")
    print("=" * 70)
    for run in runs:
        do = run['do']
        t = do.grid.t[:]
        n_re = do.eqsys.n_re.get()[:, 0]
        t_avg, gamma, gamma_norm = compute_growth_rates(do, n_re_floor)

        print(f"\n[{run['label']}]  ({run['path']})")
        print(f"  n_re initial: {n_re[0]:.3e} m^-3   final: {n_re[-1]:.3e} m^-3")

        if len(gamma) == 0:
            print("  No valid growth rate samples.")
            continue

        mask = t_avg >= t_min
        if not np.any(mask):
            mask = np.ones_like(t_avg, dtype=bool)
        g, gn, tt = gamma[mask], gamma_norm[mask], t_avg[mask]

        print(f"  Final gamma: {g[-1]:.4e} s^-1 at t = {tt[-1]:.4f} s")
        if len(g) > 10:
            avg_last = np.mean(g[int(0.9 * len(g)):])
            print(f"  Average gamma (last 10%): {avg_last:.4e} s^-1")
            lo, hi = len(g) // 10, 9 * len(g) // 10
            print(f"  Mean gamma (middle 80%): {np.mean(g[lo:hi]):.4e} s^-1"
                  f" ± {np.std(g[lo:hi]):.3e}")
        if np.any(np.isfinite(gn)):
            print(f"  Final normalized gamma: {gn[-1]:.4f}")
            if run['ef'] is not None:
                rp = rosenbluth_norm(run['ef'])
                print(f"  Rosenbluth 2*(E/Ec-1) = 2*({run['ef']:.4f}-1) = {rp:.4f}")
    print("=" * 70)


# ───────────────────────── CLI / main ─────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description='Runaway electron density and growth rate analysis')
    parser.add_argument('--data_dir', '--data-dir', type=str, dest='data_dir',
                        nargs='+', default=['../outputs/dreicer_with_fre_output.h5'],
                        help='One or more DREAM output HDF5 files. '
                             '1 file → single-run mode; several → E/Ec scan mode.')
    parser.add_argument('--plot_dir', '--plot-dir', type=str, dest='plot_dir',
                        default='../figures',
                        help='Directory to save plots')
    parser.add_argument('--t-min', type=float, default=0.5,
                        help='Start time (s) for all time-series plots (default: 0.5)')
    parser.add_argument('--n-re-floor', type=float, default=1e5,
                        help='Minimum n_re for valid growth-rate samples (default: 1e5)')
    parser.add_argument('--show', action='store_true',
                        help='Call plt.show() after saving (default: off)')
    return parser.parse_args()


def main():
    args = parse_args()

    print("=" * 70)
    print("DREAM runaway density & growth rate analysis")
    print("=" * 70)

    runs = load_runs(args.data_dir)
    if not runs:
        print("Error: no simulation results loaded!")
        return

    os.makedirs(args.plot_dir, exist_ok=True)
    single = len(runs) == 1

    if single:
        do = runs[0]['do']
        print("\nSingle-run mode: detailed density and growth rate figures")
        plot_density_evolution(do, args.plot_dir, t_min=args.t_min)
        plot_nre_growth_2panel(do, args.plot_dir, t_min=args.t_min,
                               n_re_floor=args.n_re_floor)
    else:
        print(f"\nMulti-run E/Ec scan mode: {len(runs)} runs")
        plot_scan_overlays(runs, args.plot_dir, t_min=args.t_min,
                           n_re_floor=args.n_re_floor)
        plot_avg_growth_vs_Ec(runs, args.plot_dir, n_re_floor=args.n_re_floor)

    plot_normalized_growth_rate(runs, args.plot_dir, t_min=args.t_min,
                                n_re_floor=args.n_re_floor)
    print_stats(runs, n_re_floor=args.n_re_floor, t_min=args.t_min)

    if args.show:
        plt.show()

    for run in runs:
        run['do'].close()

    print(f"\nPlots saved to {args.plot_dir}")


if __name__ == '__main__':
    main()
