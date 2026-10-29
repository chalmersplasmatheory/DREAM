#!/usr/bin/env python3
#
# Distribution function visualization for DREAM output (consolidated).
#
# Figures:
#   contour  - p-xi log(f) contour (single time, or one PNG per time in batch mode;
#              saved in f_contour/)
#   zoom     - zoomed p-xi contour with custom axis limits (saved in f_zoom/)
#   grid     - multi-time p-xi contour panels (optional, not in defaults)
#   spectra  - xi-integrated f(p) overlay at multiple times
#   detail   - 4-time detail panels of f(p)
#   pitch    - f(p) at different pitch angles (single time)
#   trapped  - trapped-electron spectrum in kinetic energy
#
# Merges plot_logf_px_contour.py + plot_f_time_evolution.py + plot_dreicer_dist.py (fig 1).
#
# Examples:
#   python plot_f.py --data_dir ../outputs/quasilinear_whistler_output.h5 --species f_re --figures all
#   python plot_f.py --t-min 1.0 --t-max 2.0 --dt 0.01 --figures contour   # one PNG per 0.01 s
#   python plot_f.py --time_series 0.9,1.5,2.0 --figures spectra,grid
#

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import DREAMOutput

# Physical constants
C = 299792458.0               # Speed of light (m/s)
E_CHARGE = 1.602176634e-19    # Elementary charge (C)
M_ELECTRON = 9.10938356e-31   # Electron mass (kg)

SPECIES_LABEL = {
    'f_re': 'Runaway electrons ($f_\\mathrm{re}$)',
    'f_hot': 'Hot electrons ($f_\\mathrm{hot}$)',
}

FIGURE_TOKENS = ('contour', 'grid', 'zoom', 'spectra', 'detail', 'pitch', 'trapped')
DEFAULT_FIGURES = ('contour', 'spectra', 'detail', 'pitch', 'trapped')


# ───────────────────────── resonance helpers ─────────────────────────

def calculate_resonance_lines(p_max, omega, k_par, B0, n_values=None):
    """
    Calculate resonance condition curves in (p, xi) space.

    Resonance equation: omega - k_par*v_par - n*Omega_ce = 0
    Solved for xi as a function of p:
        xi = (omega*gamma - n*Omega_ce) / (k_par*p*c)

    Returns {n: (p_plot, xi_plot)} with |xi| <= 1 masked arrays.
    """
    if n_values is None:
        n_values = [0, -1, 1]

    omega_ce = E_CHARGE * B0 / M_ELECTRON

    p_res = np.linspace(0.01, p_max, 2000)
    gamma_res = np.sqrt(p_res**2 + 1)

    resonance_lines = {}
    for n in n_values:
        xi = (omega * gamma_res - n * omega_ce) / (k_par * p_res * C)
        mask = np.abs(xi) <= 1
        resonance_lines[n] = (p_res[mask], xi[mask])

    return resonance_lines


def plot_resonance_on_ax(ax, resonance_lines, n_values=None):
    """Overlay resonance curves on an existing Axes."""
    if n_values is None:
        n_values = list(resonance_lines.keys())

    color_map = {
        0: 'yellow',
        1: 'magenta',
        -1: 'cyan',
        2: 'lime',
        -2: 'orange',
    }

    for n in n_values:
        if n not in resonance_lines:
            continue
        p_plot, xi_plot = resonance_lines[n]
        if len(p_plot) == 0:
            continue
        color = color_map.get(n, 'white')
        ax.plot(p_plot, xi_plot, linestyle='--', color=color, linewidth=2,
                label=f'n={n:+d} resonance')

    return ax


# ───────────────────────── data helpers ─────────────────────────

def load_species(do, species):
    """Load distribution f [nt, nxi, np] and its (p, xi) grid."""
    if species == 'f_re':
        f = do.eqsys.f_re[:, 0, :, :]
        p = do.grid.runaway.p[:]
        xi = do.grid.runaway.xi[:]
    elif species == 'f_hot':
        f = do.eqsys.f_hot[:, 0, :, :]
        p = do.grid.hottail.p[:]
        xi = do.grid.hottail.xi[:]
    else:
        raise ValueError("species must be 'f_re' or 'f_hot'")
    return f, p, xi


def has_f_re(do):
    return hasattr(do.eqsys, 'f_re') and do.eqsys.f_re is not None


def log10_f_normalized(f_data):
    """log10(|f| / f_max + 1e-30) — normalized log distribution for contours."""
    f_abs = np.abs(f_data)
    f_max = np.max(f_abs)
    if f_max > 0:
        return np.log10(f_abs / f_max + 1e-30)
    return np.zeros_like(f_abs)


def _auto_ylim(values, pad_bottom=3, pad_top=3):
    """Compute a sensible ylim from data, ignoring non-finite entries."""
    valid = np.isfinite(values)
    if not np.any(valid):
        return -20, 5
    vmin = np.min(values[valid])
    vmax = np.max(values[valid])
    return vmin - pad_bottom, vmax + pad_top


def resolve_time_indices(t_grid, time_series=None, t_min=None, t_max=None, dt=None):
    """
    Map requested times (s) to nearest stored time-step indices.

    Priority: explicit time_series > [t_min, t_max] (+ dt spacing) > 5 evenly spaced.
    Returns order-preserving de-duplicated list of int indices.
    """
    t_grid = np.asarray(t_grid)
    nt = len(t_grid)

    if time_series:
        target_times = np.asarray(time_series, dtype=float)
    elif t_min is not None and t_max is not None:
        if dt is not None and dt > 0:
            target_times = np.arange(t_min, t_max + 0.5 * dt, dt)
        else:
            target_times = np.linspace(t_min, t_max, 5)
    else:
        return [int(i) for i in np.linspace(0, nt - 1, 5, dtype=int)]

    indices = [int(np.argmin(np.abs(t_grid - tt))) for tt in target_times]
    return list(dict.fromkeys(indices))


def format_time_names(t_vals):
    """Format time values for filenames; fall back to more digits on collision."""
    names = [f'{t:.3f}' for t in t_vals]
    if len(set(names)) < len(names):
        names = [f'{t:.6f}' for t in t_vals]
    return names


def title_info(do, t_idx):
    """(a/R text, E text) annotations read from the output, or None on failure."""
    a_R_text = None
    E_text = None
    try:
        a = do.grid.a[:]
        R0 = do.grid.R0[:]
        a_val = float(a[0]) if hasattr(a, '__len__') and len(a) > 0 else float(a)
        R0_val = float(R0[0]) if hasattr(R0, '__len__') and len(R0) > 0 else float(R0)
        a_R_text = f'a/R = {a_val / R0_val:.2f}'
    except Exception as e:
        print(f"Warning: Could not read a/R value: {e}")
    try:
        E_field = do.eqsys.E_field.get()
        E_val = E_field[t_idx, 0]
        E_text = f'E = {E_val:.2f} V/m'
    except Exception as e:
        print(f"Warning: Could not read E_field: {e}")
    return a_R_text, E_text


def species_name(species):
    return SPECIES_LABEL.get(species, species)


# ───────────────────────── p-xi contour figures ─────────────────────────

def plot_pxi_contour(do, f, p, xi, t_idx, species, plot_dir, resonance_params,
                     p_range=None, xi_range=None, time_label=None):
    """
    Single p-xi log(f) contour with a/R and E annotations in the title.
    If p_range/xi_range given, produce the zoomed variant (caller should pass
    a dedicated output directory, e.g. <plot_dir>/f_zoom).
    """
    zoom = p_range is not None or xi_range is not None
    P, XI = np.meshgrid(p, xi)
    log_f = log10_f_normalized(f[t_idx, :, :])

    fig, ax = plt.subplots(1, 1, figsize=(6, 5) if zoom else (10, 8))

    levels = np.linspace(-12, 0, 25)
    contour = ax.contourf(P, XI, log_f, levels=levels, cmap='viridis', extend='both')
    ax.contour(P, XI, log_f, levels=levels[::2], colors='white',
               linewidths=0.5, alpha=0.5)

    if resonance_params is not None:
        p_max = p_range[1] if (zoom and p_range is not None) else np.max(p)
        rl = calculate_resonance_lines(
            p_max,
            omega=resonance_params['omega'],
            k_par=resonance_params['k_par'],
            B0=resonance_params['B0'],
            n_values=resonance_params.get('n_values', None),
        )
        plot_resonance_on_ax(ax, rl, n_values=resonance_params.get('n_values', None))

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)

    ax.set_xlabel(r'Momentum $p$ [$m_e c$]', fontsize=12)
    ax.set_ylabel(r'Pitch $\xi = \cos(\theta)$', fontsize=12)

    if zoom:
        if p_range is not None:
            ax.set_xlim(p_range[0], p_range[1])
        if xi_range is not None:
            ax.set_ylim(xi_range[0], xi_range[1])
    else:
        ax.invert_yaxis()

    t_val = do.grid.t[t_idx]
    a_R_text, E_text = title_info(do, t_idx)
    title = f'{species_name(species)}'
    if zoom:
        title += ' (zoomed)'
    title += f'\nTime = {t_val:.6e} s'
    if a_R_text:
        title += f'\n{a_R_text}'
    if E_text:
        title += f', {E_text}'
    ax.set_title(title, fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.3)
    if resonance_params is not None:
        ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    os.makedirs(plot_dir, exist_ok=True)
    if time_label is None:
        time_label = f'{t_val:.3f}'
    fname = f'{species}_t{time_label}.png'
    save_path = os.path.join(plot_dir, fname)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close(fig)
    return save_path


def plot_pxi_grid(do, f, p, xi, t_indices, species, plot_dir, resonance_params):
    """Multi-time p-xi contour panels sharing one colorbar."""
    P, XI = np.meshgrid(p, xi)

    resonance_lines = None
    if resonance_params is not None:
        resonance_lines = calculate_resonance_lines(
            np.max(p),
            omega=resonance_params['omega'],
            k_par=resonance_params['k_par'],
            B0=resonance_params['B0'],
            n_values=resonance_params.get('n_values', None),
        )

    n_plots = len(t_indices)
    fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 5))
    if n_plots == 1:
        axes = [axes]
    else:
        axes = list(axes)

    contour = None
    for idx, t_idx in enumerate(t_indices):
        log_f = log10_f_normalized(f[t_idx, :, :])
        levels = np.linspace(-12, 0, 25)
        contour = axes[idx].contourf(P, XI, log_f, levels=levels,
                                     cmap='viridis', extend='both')

        axes[idx].set_xlabel(r'$p$ [$m_e c$]', fontsize=10)
        if idx == 0:
            axes[idx].set_ylabel(r'$\xi$', fontsize=10)
        axes[idx].invert_yaxis()

        t_val = do.grid.t[t_idx]
        axes[idx].set_title(f't = {t_val:.6e} s', fontsize=11)
        axes[idx].grid(True, linestyle='--', alpha=0.3)

        if resonance_lines is not None:
            plot_resonance_on_ax(axes[idx], resonance_lines,
                                 n_values=resonance_params.get('n_values', None))
            axes[idx].legend(loc='upper right', fontsize=8)

    cbar = fig.colorbar(contour, ax=axes, fraction=0.02, pad=0.04)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)

    fig.suptitle(f'{species_name(species)} distribution evolution',
                 fontsize=14, y=1.02)
    plt.tight_layout()

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, f'logf_evolution_{species}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close(fig)
    return save_path


# ───────────────────────── 1D spectrum figures ─────────────────────────

def plot_spectra_overlay(do, f, p, xi, t_indices, species, plot_dir):
    """xi-integrated f(p) curves overlaid at multiple times."""
    plt.figure(figsize=(6, 6))

    line_colors = ['#1f77b4', '#ff6b6b', '#2ca02c', '#ff7f0e', '#9467bd']
    line_styles = ['--'] + ['-'] * (len(t_indices) - 1)
    all_log = []

    for i, t_idx in enumerate(t_indices):
        f_at_time = f[t_idx, :, :]
        f_vs_p = np.trapz(np.abs(f_at_time), xi, axis=0)
        f_nonzero = np.where(f_vs_p > 0, f_vs_p, 1e-30)
        log_f = np.log10(f_nonzero)
        all_log.append(log_f)

        plt.plot(p, log_f, color=line_colors[i % len(line_colors)],
                 linestyle=line_styles[i % len(line_styles)], linewidth=2,
                 label=f't={do.grid.t[t_idx]:.3f} s')

    plt.xlabel('p ($m_e c$)')
    plt.ylabel(r'$\log_{10}( \int |f| \, d\xi )$')
    plt.xlim(0, 20)
    yl, yh = _auto_ylim(np.concatenate(all_log))
    plt.ylim(yl, yh)
    plt.title(rf'{species_name(species)}: $\xi$-integrated distribution evolution')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, f'f_time_evolution_{species}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close()
    return save_path


def plot_spectra_detail(do, f, p, xi, species, plot_dir):
    """4 evenly spaced time slices of the xi-integrated f(p) (skip t=0)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    num_times = f.shape[0]
    detail_times = np.linspace(1, num_times - 1, 4, dtype=int)

    for i, t_idx in enumerate(detail_times):
        ax = axes[i]
        f_at_time = f[t_idx, :, :]
        f_vs_p = np.trapz(np.abs(f_at_time), xi, axis=0)
        f_nonzero = np.where(f_vs_p > 0, f_vs_p, 1e-30)
        log_f = np.log10(f_nonzero)

        ax.plot(p, log_f, 'b-', linewidth=2)
        ax.set_xlabel('p ($m_e c$)')
        ax.set_ylabel(r'$\log_{10}(\int |f| \, d\xi)$')
        ax.set_xlim(0, 50)
        yl, yh = _auto_ylim(log_f)
        ax.set_ylim(yl, yh)
        ax.set_title(f't = {do.grid.t[t_idx]*1e6:.1f} μs')
        ax.grid(True, linestyle='--', alpha=0.7)

    fig.suptitle(f'{species_name(species)}: distribution at 4 times', fontsize=14)
    plt.tight_layout()

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, f'f_time_detail_{species}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close(fig)
    return save_path


def plot_pitch_comparison(do, f, p, xi, t_idx, species, plot_dir):
    """f(p) at several pitch angles at one time."""
    plt.figure(figsize=(6, 6))
    nxi = xi.size
    xi_indices = np.linspace(0, nxi - 1, 5, dtype=int)
    colors = plt.cm.plasma(np.linspace(0, 1, len(xi_indices)))

    all_log = []
    for i, idx in enumerate(xi_indices):
        f_at_xi = f[t_idx, idx, :]
        f_abs = np.abs(f_at_xi)
        f_nonzero = np.where(f_abs > 0, f_abs, 1e-30)
        log_f = np.log10(f_nonzero)
        all_log.append(log_f)

        plt.plot(p, log_f, color=colors[i], linewidth=2,
                 label=f'ξ={xi[idx]:.2f}')

    plt.xlabel('p ($m_e c$)')
    plt.ylabel(r'$\log_{10}(|f|)$')
    plt.xlim(0, 50)
    yl, yh = _auto_ylim(np.concatenate(all_log))
    plt.ylim(yl, yh)
    plt.title(rf'{species_name(species)}: distribution at different pitch angles '
              f'(t={do.grid.t[t_idx]:.3f} s)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, f'f_pitch_comparison_{species}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close()
    return save_path


def plot_trapped_spectrum(do, f, p, xi, t_indices, species, plot_dir):
    """Trapped-electron (|xi| < xi0) spectrum in kinetic energy, multiple times."""
    if not (hasattr(do.grid, 'xi0TrappedBoundary') and do.grid.xi0TrappedBoundary is not None):
        print("Warning: xi0TrappedBoundary not available; skipping trapped spectrum.")
        return None

    xi0 = do.grid.xi0TrappedBoundary[0]
    mask_trapped = np.abs(xi) < xi0

    # Kinetic energy [MeV] = (sqrt(1+p^2) - 1) * 0.511
    energy = (np.sqrt(1 + p**2) - 1) * 0.511

    plt.figure(figsize=(6, 6))
    colors = plt.cm.turbo(np.linspace(0, 1, len(t_indices)))
    all_log = []

    for i, t_idx in enumerate(t_indices):
        f_at_time = f[t_idx, :, :]
        f_trapped = np.trapz(np.abs(f_at_time[mask_trapped, :]), xi[mask_trapped], axis=0)
        f_nonzero = np.where(f_trapped > 0, f_trapped, 1e-30)
        log_f = np.log10(f_nonzero)
        all_log.append(log_f)

        plt.plot(energy, log_f, color=colors[i], linewidth=2,
                 label=f't={do.grid.t[t_idx]:.3f} s')

    plt.xlabel(r'Kinetic energy $E_k$ (MeV)')
    plt.ylabel(r'$\log_{10}( \int_{\rm trapped} |f| \, d\xi )$')
    plt.xlim(0, None)
    yl, yh = _auto_ylim(np.concatenate(all_log))
    plt.ylim(yl, yh)
    plt.title(rf'{species_name(species)}: trapped-electron distribution '
              r'($|\xi| < \xi_0$)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    os.makedirs(plot_dir, exist_ok=True)
    save_path = os.path.join(plot_dir, f'f_trapped_energy_{species}.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    plt.close()
    return save_path


# ───────────────────────── CLI / main ─────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot DREAM distribution functions (p-xi contours and 1D spectra)')
    parser.add_argument('--data_dir', '--data-dir', type=str, dest='data_dir',
                        default='../outputs/dreicer_with_fre_output.h5',
                        help='Path to the DREAM output HDF5 file')
    parser.add_argument('--plot_dir', '--plot-dir', type=str, dest='plot_dir',
                        default='../figures',
                        help='Directory to save plots')
    parser.add_argument('--species', type=str, default='f_re', choices=['f_re', 'f_hot'],
                        help='Which distribution to plot (default: f_re)')
    parser.add_argument('--t', type=float, default=None,
                        help='Time (s) for single-time figures (contour/zoom/pitch); '
                             'default: last time step')
    parser.add_argument('--t-min', type=float, default=None,
                        help='Start time (s) for batch contour mode / multi-time selection')
    parser.add_argument('--t-max', type=float, default=None,
                        help='End time (s) for batch contour mode / multi-time selection')
    parser.add_argument('--dt', type=float, default=None,
                        help='Time interval (s) between batch contour plots, e.g. 0.01')
    parser.add_argument('--figures', type=str,
                        default=','.join(DEFAULT_FIGURES),
                        help='Comma-separated figure list: contour,grid,zoom,spectra,detail,'
                             'pitch,trapped,all (default: contour,grid,spectra,detail,pitch,trapped)')
    parser.add_argument('--time_series', '--time-series', type=str, dest='time_series',
                        default=None,
                        help='Comma-separated times (s) for multi-time figures, '
                             'e.g. "0.5,0.9,1.1,1.2,3.0"')
    parser.add_argument('--omega', type=float, default=2 * np.pi * 476e6,
                        help='Wave angular frequency (rad/s), default 2π×476 MHz')
    parser.add_argument('--k-par', '--k_par', type=float, dest='k_par', default=-41.0,
                        help='Parallel wavenumber (m^-1), default -41')
    parser.add_argument('--B0', type=float, default=1.4,
                        help='Magnetic field (T), default 1.4')
    parser.add_argument('--harmonics', type=str, default='-1,0,1,2',
                        help='Comma-separated resonance harmonic numbers (default: -1,0,1,2)')
    parser.add_argument('--no-resonance', action='store_true',
                        help='Disable resonance curve overlay')
    parser.add_argument('--xi-max', type=float, default=1.0,
                        help='Maximum xi for zoomed view (default: 1.0)')
    parser.add_argument('--xi-min', type=float, default=0.25,
                        help='Minimum xi for zoomed view (default: 0.25)')
    parser.add_argument('--p-max', type=float, default=30,
                        help='Maximum p for zoomed view (default: 30)')
    parser.add_argument('--p-min', type=float, default=0,
                        help='Minimum p for zoomed view (default: 0)')
    parser.add_argument('--show', action='store_true',
                        help='Call plt.show() after saving (default: off)')
    return parser.parse_args()


def main():
    args = parse_args()

    figure_tokens = [tok.strip() for tok in args.figures.split(',') if tok.strip()]
    if 'all' in figure_tokens:
        figure_tokens = list(DEFAULT_FIGURES) + ['zoom']
    unknown = [tok for tok in figure_tokens if tok not in FIGURE_TOKENS]
    if unknown:
        print(f"Error: unknown figure tokens: {unknown}. Valid: {FIGURE_TOKENS}, all")
        return
    figures = set(figure_tokens)

    if not os.path.exists(args.data_dir):
        print(f"Error: output file {args.data_dir} not found!")
        return

    print(f"Loading data from {args.data_dir}...")
    do = DREAMOutput(args.data_dir)

    species = args.species
    if species == 'f_re' and not has_f_re(do):
        print("Warning: f_re not found in output; falling back to f_hot.")
        species = 'f_hot'
    f, p, xi = load_species(do, species)
    t_grid = do.grid.t[:]
    nt = len(t_grid)

    # Resonance parameters
    resonance_params = None
    if not args.no_resonance:
        n_values = [int(x.strip()) for x in args.harmonics.split(',')]
        resonance_params = {
            'omega': args.omega,
            'k_par': args.k_par,
            'B0': args.B0,
            'n_values': n_values,
        }
        print(f"  Resonance overlay: omega={args.omega:.3e} rad/s, "
              f"k_par={args.k_par} m^-1, B0={args.B0} T, n={n_values}")
    else:
        print("  Resonance overlay: disabled")

    # Time selections
    time_series = None
    if args.time_series:
        time_series = [float(x.strip()) for x in args.time_series.split(',')]

    multi_indices = resolve_time_indices(t_grid, time_series, args.t_min, args.t_max, args.dt)

    if args.t is not None:
        single_idx = int(np.argmin(np.abs(t_grid - args.t)))
    else:
        single_idx = nt - 1

    batch_mode = args.t_min is not None and args.t_max is not None
    if batch_mode:
        # Batch contour spacing comes from --t-min/--t-max/--dt only,
        # independent of --time_series (which drives the multi-time figures).
        contour_indices = resolve_time_indices(t_grid, None, args.t_min, args.t_max, args.dt)
    else:
        contour_indices = [single_idx]
    contour_labels = format_time_names([t_grid[i] for i in contour_indices])

    print(f"  Species: {species}")
    print(f"  Single-time figures at t = {t_grid[single_idx]:.6e} s (index {single_idx})")
    print(f"  Multi-time figures at indices: {multi_indices}")
    if batch_mode:
        print(f"  Batch contour mode: {len(contour_indices)} figures "
              f"(t = [{args.t_min}, {args.t_max}], dt = {args.dt})")

    os.makedirs(args.plot_dir, exist_ok=True)

    # ── p-xi contours ──
    if 'contour' in figures:
        # Batch contour series go to a separate f_contour/ subfolder
        contour_dir = os.path.join(args.plot_dir, 'f_contour')
        for t_idx, label in zip(contour_indices, contour_labels):
            plot_pxi_contour(do, f, p, xi, t_idx, species, contour_dir,
                             resonance_params, time_label=label)
            if args.show:
                plt.show()

    if 'zoom' in figures:
        # Zoomed views go to a separate f_zoom/ subfolder
        zoom_dir = os.path.join(args.plot_dir, 'f_zoom')
        p_range = (args.p_min, args.p_max)
        xi_range = (args.xi_max, args.xi_min)
        for t_idx, label in zip(contour_indices, contour_labels):
            plot_pxi_contour(do, f, p, xi, t_idx, species, zoom_dir,
                             resonance_params, p_range=p_range, xi_range=xi_range,
                             time_label=label)
            if args.show:
                plt.show()

    if 'grid' in figures:
        plot_pxi_grid(do, f, p, xi, multi_indices, species, args.plot_dir,
                      resonance_params)
        if args.show:
            plt.show()

    # ── 1D spectra ──
    if 'spectra' in figures:
        plot_spectra_overlay(do, f, p, xi, multi_indices, species, args.plot_dir)
        if args.show:
            plt.show()

    if 'detail' in figures:
        plot_spectra_detail(do, f, p, xi, species, args.plot_dir)
        if args.show:
            plt.show()

    if 'pitch' in figures:
        plot_pitch_comparison(do, f, p, xi, single_idx, species, args.plot_dir)
        if args.show:
            plt.show()

    if 'trapped' in figures:
        plot_trapped_spectrum(do, f, p, xi, multi_indices, species, args.plot_dir)
        if args.show:
            plt.show()

    print(f"Plots saved to {args.plot_dir}")
    do.close()


if __name__ == '__main__':
    main()
