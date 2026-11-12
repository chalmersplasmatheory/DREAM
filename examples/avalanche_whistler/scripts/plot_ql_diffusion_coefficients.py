"""
Visualize DREAM quasilinear diffusion coefficients in phase space
and compute the injected wave energy density.

Diffusion coefficients
----------------------
The real net FVM diffusion coefficients on the runaway grid
('runaway/Dpp', 'runaway/Dpx', 'runaway/Dxp', 'runaway/Dxx') are read from
the output file (enable them via ds.other.include in the generator).
They are shown as a 4-panel phase-space snapshot at the time of maximum
|D_pp| (i.e. during wave injection), together with the time evolution of
the diffusion strength. If the quantities are missing from the output,
no figure is produced (no fake/demo data).

Wave energy
-----------
The normalized wave energy density per mode is the 'amplitude' parameter of
the quasilinear settings:

    E_k = |E_kx|^2 * eps0 * D / (4 m_e c^2)      ( = amplitude )

so the total wave energy density is

    W = m_e c^2 \\int d^3k E_k
      = 2 pi m_e c^2 \\iint E_k k^2 sin(theta_k) dk dtheta_k
      = 2 pi m_e c^2 E_k (k_max^3 - k_min^3)/3 (cos(theta_min) - cos(theta_max))

for a uniform spectrum. The time dependence W(t) follows the injection
envelope A(t) (start time, 50% duty cycle, ramp-up).

Usage:
    Run in the scripts directory: python plot_ql_diffusion_coefficients.py [options]

Command-line arguments:
    --data_dir      Data directory (default: ../outputs)
    --plot_dir      Plot output directory (default: ../figures)
    --output_file   Output file name (default: quasilinear_whistler_output.h5)
    --t             Time (s) for the phase-space snapshot
                    (default: the time of maximum |D_pp|)

Examples:
    python plot_ql_diffusion_coefficients.py
    python plot_ql_diffusion_coefficients.py --data_dir ../outputs --plot_dir ../figures_w
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py
import os
import sys
import argparse

# Add DREAM Python library path
sys.path.insert(0, '/data/zhzhou/DREAM/py')
from DREAM import DREAMOutput

# Physical constants
M_E = 9.10938356e-31        # Electron mass (kg)
C_LIGHT = 299792458.0       # Speed of light (m/s)


# ───────────────────────── wave energy ─────────────────────────

def load_ql_wave_settings(filepath):
    """
    Read the quasilinear wave settings (amplitude, spectrum, injection timing)
    from a DREAM output file (settings/eqsys/...) or settings file (eqsys/...).

    Returns a dict, or None if no quasilinear settings are present.
    """
    if not os.path.exists(filepath):
        print(f"Error: file not found {filepath}")
        return None

    candidates = [
        'settings/eqsys/f_re/quasilinear',
        'settings/eqsys/f_hot/quasilinear',
        'eqsys/f_re/quasilinear',
        'eqsys/f_hot/quasilinear',
    ]
    with h5py.File(filepath, 'r') as f:
        for path in candidates:
            if path in f:
                g = f[path]
                out = {}
                for key in g.keys():
                    val = g[key][()]
                    out[key] = val.item() if hasattr(val, 'item') and val.size == 1 else val
                print(f"Quasilinear wave settings found at '{path}'")
                return out
    print("Warning: no quasilinear settings found in the file")
    return None


def wave_energy_density(eps_k, k_min, k_max, ktheta_min, ktheta_max):
    """
    Total wave energy density from the normalized per-mode energy density
    E_k (= DREAM 'amplitude') for a uniform spectrum over
    [k_min, k_max] x [ktheta_min, ktheta_max]:

        W = m_e c^2 \\int d^3k E_k
          = 2 pi m_e c^2 E_k (k_max^3 - k_min^3)/3 (cos ktheta_min - cos ktheta_max)

    Returns W in J/m^3 (SI) for k in m^-1.
    """
    k_space_volume = (2.0 * np.pi
                      * (k_max**3 - k_min**3) / 3.0
                      * (np.cos(ktheta_min) - np.cos(ktheta_max)))
    return eps_k * M_E * C_LIGHT**2 * k_space_volume


def injection_envelope(t, amplitude, start_inject_time, inject_cycle_duration, ramp_time):
    """
    Effective wave amplitude A(t) of the periodic injection schedule
    (mirrors QuasilinearDiffusionTerm::calculateEffectiveAmplitude):
    - t < start_inject_time: A = 0
    - inject_cycle_duration <= 0: continuous injection
    - else: 50% duty-cycle ON/OFF, linear ramp-up over ramp_time at each ON start
    """
    t = np.asarray(t, dtype=float)
    if start_inject_time is not None and start_inject_time >= 0:
        active = t >= start_inject_time
        t_local = t - start_inject_time
    else:
        active = np.ones_like(t, dtype=bool)
        t_local = t

    if not inject_cycle_duration or inject_cycle_duration <= 0:
        on = active
        phase = t_local
    else:
        phase = np.mod(t_local, inject_cycle_duration)
        on = active & (phase < inject_cycle_duration / 2.0)

    if ramp_time and ramp_time > 0:
        ramp = np.minimum(phase / ramp_time, 1.0)
    else:
        ramp = np.ones_like(t)

    return np.where(on, amplitude * ramp, 0.0)


def compute_injected_wave_energy(filepath, plot_dir, output_prefix='ql_diffusion'):
    """
    Compute the injected wave energy density W from the quasilinear settings:
    peak value (full amplitude) and W(t) along the injection envelope.
    Prints a summary and saves a CSV + plot of W(t).
    """
    qs = load_ql_wave_settings(filepath)
    if qs is None:
        return None

    amplitude = float(qs.get('amplitude', 0.0))
    k_min = float(qs.get('k_min', 0.0))
    k_max = float(qs.get('k_max', 0.0))
    ktheta_min = float(qs.get('ktheta_min', 0.0))
    ktheta_max = float(qs.get('ktheta_max', 0.0))
    num_k = int(qs.get('num_k', 0))
    num_ktheta = int(qs.get('num_ktheta', 0))
    start_inject_time = float(qs.get('start_inject_time', -1.0))
    inject_cycle_duration = float(qs.get('inject_cycle_duration', 0.0))
    ramp_time = float(qs.get('ramp_time', 0.0))

    W_peak = wave_energy_density(amplitude, k_min, k_max, ktheta_min, ktheta_max)

    print("\n" + "=" * 70)
    print("INJECTED WAVE ENERGY DENSITY")
    print("=" * 70)
    print("  Normalized wave energy density per mode: E_k = amplitude")
    print(f"    amplitude      = {amplitude:.4e}")
    print(f"    spectrum       = uniform {num_k} x {num_ktheta} modes")
    print(f"    k range        = [{k_min:.3f}, {k_max:.3f}] m^-1")
    print(f"    theta_k range  = [{ktheta_min:.3f}, {ktheta_max:.3f}] rad")
    k_space_volume = (2.0 * np.pi
                      * (k_max**3 - k_min**3) / 3.0
                      * (np.cos(ktheta_min) - np.cos(ktheta_max)))
    print(f"    k-space factor = 2*pi*(k_max^3-k_min^3)/3*(cos th_min - cos th_max)"
          f" = {k_space_volume:.4e} m^-3")
    print(f"  Peak wave energy density (full amplitude):")
    print(f"    W_peak = m_e c^2 * E_k * k-space factor = {W_peak:.6e} J/m^3")
    print(f"  Injection: start={start_inject_time} s, cycle={inject_cycle_duration} s, "
          f"ramp={ramp_time} s")

    # Time dependence along the injection envelope
    with h5py.File(filepath, 'r') as f:
        t = None
        if 'grid/t' in f:
            t = f['grid/t'][:]
        elif 'timestep/tmax' in f:
            t_max = float(f['timestep/tmax'][()])
            t = np.linspace(0.0, t_max, 2000)

    if t is not None:
        A_t = injection_envelope(t, amplitude, start_inject_time,
                                 inject_cycle_duration, ramp_time)
        W_t = wave_energy_density(A_t, k_min, k_max, ktheta_min, ktheta_max)

        os.makedirs(plot_dir, exist_ok=True)

        csv_path = os.path.join(plot_dir, f'{output_prefix}_wave_energy.csv')
        with open(csv_path, 'w') as fh:
            fh.write("time,amplitude,W_J_per_m3\n")
            for ti, ai, wi in zip(t, A_t, W_t):
                fh.write(f"{ti:.6e},{ai:.6e},{wi:.6e}\n")
        print(f"  Wave energy time series saved to: {csv_path}")

        plt.figure(figsize=(7, 5))
        plt.plot(t, W_t, 'b-', linewidth=2)
        plt.xlabel('Time (s)')
        plt.ylabel(r'Wave energy density $\overline{W}$ (J/m$^3$)')
        plt.title(f'Injected wave energy density vs time\n'
                  f'$W_{{peak}}$ = {W_peak:.3e} J/m$^3$, '
                  f'start={start_inject_time} s, cycle={inject_cycle_duration} s')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plot_path = os.path.join(plot_dir, f'{output_prefix}_wave_energy_vs_time.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Figure saved: {plot_path}")
    else:
        print("  Warning: no time grid found; skipping W(t) output")

    print("=" * 70)
    return W_peak


# ───────────────────────── diffusion coefficients ─────────────────────────

def load_diffusion_coefficients(filepath, grid='runaway'):
    """
    Read the net FVM diffusion coefficients ('<grid>/Dpp', 'Dpx', 'Dxp', 'Dxx')
    from a DREAM output file and average them onto momentum-grid cell centers.

    Returns dict with 't', 'p', 'xi' and 3D arrays (nt, np1, np2) for each of
    Dpp/Dpx/Dxp/Dxx, or None if the quantities are not stored.
    """
    quants = {}
    with h5py.File(filepath, 'r') as f:
        for short in ('Dpp', 'Dpx', 'Dxp', 'Dxx'):
            path = f'{grid}/{short}'
            if f'other/{path}' not in f:
                print(f"Warning: 'other/{path}' not found in {filepath}")
                return None
            quants[short] = f[f'other/{path}'][:]

        p = f[f'grid/{grid}/p1'][:]
        p_f = f[f'grid/{grid}/p1_f'][:]
        xi = f[f'grid/{grid}/p2'][:]
        xi_f = f[f'grid/{grid}/p2_f'][:]
        t = f['grid/t'][:]

    np1, np2 = len(p), len(xi)
    nt = quants['Dpp'].shape[0]

    # Other quantities are stored at save steps only, i.e. corresponding to
    # t[1:] when the initial state is also in grid/t.
    if len(t) == nt + 1:
        t = t[1:]
    elif len(t) != nt:
        print(f"Warning: {nt} diffusion frames vs {len(t)} time points; "
              f"aligning to the last {nt} times")
        t = t[-nt:]

    out = {'t': t, 'p': p, 'xi': xi}

    # Stored layout per frame (C order):
    #   Dpp, Dpx: (nxi, np1+1)  on p-faces
    #   Dxp, Dxx: (nxi+1, np1)  on xi-faces
    # Average onto cell centers -> (nxi, np1), ready for pcolormesh(p, xi, C).
    def to_centers_p(arr):      # (nt, nr, np2, np1+1) -> (nt, np2, np1)
        raw = arr[:, 0]
        return 0.5 * (raw[..., :-1] + raw[..., 1:])

    def to_centers_x(arr):      # (nt, nr, np2+1, np1) -> (nt, np2, np1)
        raw = arr[:, 0]
        return 0.5 * (raw[:, :-1, :] + raw[:, 1:, :])

    out['Dpp'] = to_centers_p(quants['Dpp'])
    out['Dpx'] = to_centers_p(quants['Dpx'])
    out['Dxp'] = to_centers_x(quants['Dxp'])
    out['Dxx'] = to_centers_x(quants['Dxx'])
    return out


def plot_diffusion_snapshot(data, t_idx, plot_dir, output_prefix='ql_diffusion'):
    """4-panel phase-space snapshot of the real diffusion coefficients."""
    p, xi = data['p'], data['xi']
    t_val = data['t'][t_idx]
    Dpp = data['Dpp'][t_idx]
    Dpx = data['Dpx'][t_idx]
    Dxx = data['Dxx'][t_idx]

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Dxx/Dpx span ~20 decades and are essentially zero where no resonance
    # exists; clamping the log colour scale at 1e-5 keeps the plot readable
    # (values below the floor are drawn with the lowest colour).
    panels = [
        (axes[0, 0], np.abs(Dpp), 'viridis', r'$D_{pp}$ Distribution', r'$D_{pp}$', None),
        (axes[0, 1], np.abs(Dxx), 'plasma', r'$D_{\xi\xi}$ Distribution', r'$D_{\xi\xi}$', 1e-5),
        (axes[1, 0], np.abs(Dpx), 'inferno', r'$D_{p\xi}$ Distribution', r'$D_{p\xi}$', 1e-5),
    ]
    for ax, D_plot, cmap, title, label, vmin in panels:
        D_plot = np.where(D_plot > 0, D_plot, 1e-30)
        vmax = np.max(D_plot)
        if vmin is None:
            vmin = np.min(D_plot)
        else:
            # Guard against an all-zero / pre-injection snapshot
            vmax = max(vmax, vmin * 10.0)
        im = ax.pcolormesh(p, xi, D_plot, cmap=cmap,
                           norm=LogNorm(vmin=vmin, vmax=vmax),
                           shading='nearest')
        ax.set_xlabel('p (momentum)', fontsize=12)
        ax.set_ylabel(r'$\xi$ (pitch angle cosine)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        plt.colorbar(im, ax=ax, label=label)
        ax.grid(True, alpha=0.3)

    # Correlation ratio D_pxi / sqrt(D_pp * D_xixi), in [-1, 1] (Cauchy-Schwarz)
    ax4 = axes[1, 1]
    denominator = np.sqrt(np.abs(Dpp * Dxx))
    denominator = np.where(denominator > 0, denominator, 1e-30)
    ratio = Dpx / denominator
    im4 = ax4.pcolormesh(p, xi, ratio, cmap='RdBu_r', vmin=-1, vmax=1, shading='nearest')
    ax4.set_xlabel('p (momentum)', fontsize=12)
    ax4.set_ylabel(r'$\xi$ (pitch angle cosine)', fontsize=12)
    ax4.set_title(r'$D_{p\xi} / \sqrt{D_{pp} D_{\xi\xi}}$', fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im4, ax=ax4)
    cbar.set_label('Correlation Coefficient')
    ax4.grid(True, alpha=0.3)

    # Flip the pitch-angle axis on every panel: xi = -1 at the top, +1 at the
    # bottom (more intuitive for the trapped/passing structure of f(p, xi)).
    for ax in axes.flat:
        ax.invert_yaxis()

    fig.suptitle(f'Net diffusion coefficients (runaway grid) at t = {t_val:.4f} s '
                 f'(during wave injection)', fontsize=15)
    plt.tight_layout()

    os.makedirs(plot_dir, exist_ok=True)
    output_path = os.path.join(plot_dir, f'{output_prefix}_coefficients.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Figure saved: {output_path}")
    return output_path


def plot_diffusion_strength_vs_time(data, plot_dir, qs=None,
                                    output_prefix='ql_diffusion'):
    """Time evolution of max|D_pp| and max|D_xx| (diffusion strength)."""
    t = data['t']
    dpp_max = np.max(np.abs(data['Dpp']), axis=(1, 2))
    dxx_max = np.max(np.abs(data['Dxx']), axis=(1, 2))

    plt.figure(figsize=(7, 5))
    plt.semilogy(t, np.where(dpp_max > 0, dpp_max, np.nan), 'b-', linewidth=2,
                 label=r'max $|D_{pp}|$')
    plt.semilogy(t, np.where(dxx_max > 0, dxx_max, np.nan), 'r-', linewidth=2,
                 label=r'max $|D_{\xi\xi}|$')

    if qs is not None and float(qs.get('start_inject_time', -1)) >= 0:
        t0 = float(qs['start_inject_time'])
        plt.axvline(t0, color='k', linestyle='--', linewidth=1,
                    label=f'wave injection start ({t0} s)')

    plt.xlabel('Time (s)')
    plt.ylabel('Diffusion coefficient (SI)')
    plt.title('Diffusion strength vs time\n'
              r'max $|D_{pp}|$, max $|D_{\xi\xi}|$ on the runaway grid')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    os.makedirs(plot_dir, exist_ok=True)
    output_path = os.path.join(plot_dir, f'{output_prefix}_strength_vs_time.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Figure saved: {output_path}")
    return output_path


# ───────────────────────── CLI / main ─────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Visualize DREAM quasilinear diffusion coefficients and wave energy')
    parser.add_argument('--data_dir', type=str, default='../outputs', help='Data directory')
    parser.add_argument('--plot_dir', type=str, default='../figures', help='Plot output directory')
    parser.add_argument('--output_file', type=str, default='quasilinear_whistler_output.h5',
                        help='Output file name')
    parser.add_argument('--output_prefix', type=str, default='ql_diffusion',
                        help='Output file prefix')
    parser.add_argument('--t', type=float, default=None,
                        help='Time (s) for the phase-space snapshot '
                             '(default: time of maximum |D_pp|)')

    args = parser.parse_args()

    print("="*70)
    print("DREAM quasilinear diffusion coefficient visualization")
    print("="*70)

    filepath = os.path.join(args.data_dir, args.output_file)

    # Wave energy density from the quasilinear settings
    qs = load_ql_wave_settings(filepath)
    if qs is not None:
        compute_injected_wave_energy(filepath, args.plot_dir, args.output_prefix)

    # Real diffusion coefficients from the output
    data = load_diffusion_coefficients(filepath, grid='runaway')
    if data is None:
        print("\nNo diffusion coefficients in this output.")
        print("The quantities 'runaway/Dpp', 'runaway/Dpx', 'runaway/Dxp', 'runaway/Dxx'")
        print("must be enabled in the generator (ds.other.include) and the simulation")
        print("re-run. The current generator versions already do this.")
        return

    # Pick snapshot time: requested time, or the time of maximum |D_pp|
    if args.t is not None:
        t_idx = int(np.argmin(np.abs(data['t'] - args.t)))
    else:
        dpp_max_t = np.max(np.abs(data['Dpp']), axis=(1, 2))
        t_idx = int(np.argmax(dpp_max_t))
    print(f"\nSnapshot at t = {data['t'][t_idx]:.4f} s (index {t_idx})")
    print(f"  max|Dpp| = {np.max(np.abs(data['Dpp'][t_idx])):.4e}")
    print(f"  max|Dxx| = {np.max(np.abs(data['Dxx'][t_idx])):.4e}")

    plot_diffusion_snapshot(data, t_idx, args.plot_dir, args.output_prefix)
    plot_diffusion_strength_vs_time(data, args.plot_dir, qs, args.output_prefix)

    print("\nDone!")


if __name__ == '__main__':
    main()
