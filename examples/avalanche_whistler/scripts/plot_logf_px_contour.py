#!/usr/bin/env python3
#
# Plot log(f) contour in p-xi space for runaway electrons (f_re) and hot electrons (f_hot)
#

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import DREAMOutput

# Physical constants
C = 299792458.0        # Speed of light (m/s)
E_CHARGE = 1.602176634e-19  # Elementary charge (C)
M_ELECTRON = 9.10938356e-31 # Electron mass (kg)


def calculate_resonance_lines(p_max, omega, k_par, B0, n_values=None):
    """
    Calculate resonance condition curves in (p, xi) space.

    Resonance equation: ω - k_∥·v_∥ - n·Ω_ce = 0
    Solved for xi as a function of p:
        xi = (ω·γ - n·Ω_ce) / (k_∥·p·c)

    Parameters
    ----------
    p_max : float
        Maximum momentum (in m_e*c units) for the calculation grid.
    omega : float
        Wave angular frequency (rad/s).
    k_par : float
        Parallel wavenumber (m^-1).
    B0 : float
        Magnetic field strength (T).
    n_values : list of int, optional
        Resonance harmonic numbers. Default: [0, -1, +1].

    Returns
    -------
    dict
        {n: (p_plot, xi_plot)} where p_plot and xi_plot are masked arrays
        satisfying |xi| <= 1.
    """
    if n_values is None:
        n_values = [0, -1, 1]

    # Electron cyclotron frequency (rad/s)
    omega_ce = E_CHARGE * B0 / M_ELECTRON

    # Momentum grid for resonance calculation (dense, to get smooth curves)
    p_res = np.linspace(0.01, p_max, 2000)
    gamma_res = np.sqrt(p_res**2 + 1)

    resonance_lines = {}
    for n in n_values:
        # xi = (ω·γ - n·Ω_ce) / (k_∥·p·c)
        xi = (omega * gamma_res - n * omega_ce) / (k_par * p_res * C)

        # Physical constraint: |xi| <= 1
        mask = np.abs(xi) <= 1
        resonance_lines[n] = (p_res[mask], xi[mask])

    return resonance_lines


def plot_resonance_on_ax(ax, resonance_lines, n_values=None):
    """
    Overlay resonance curves on an existing Axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw on.
    resonance_lines : dict
        Output from calculate_resonance_lines().
    n_values : list of int, optional
        Which harmonics to plot. Default: all keys in resonance_lines.
    """
    if n_values is None:
        n_values = list(resonance_lines.keys())

    # Colour mapping for different harmonic orders
    color_map = {
        0:  'yellow',
        1:  'magenta',
        -1: 'cyan',
        2:  'lime',
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


def plot_logf_contour(output_file='../outputs/dreicer_with_fre_output.h5', 
                      time_index=-1, species='f_re', plot_dir='../figures',
                      resonance_params=None):
    """
    Plot log(f) contour in p-xi space
    
    Parameters:
    -----------
    output_file : str
        Path to the DREAM output file
    time_index : int
        Time index to plot (-1 for last time step)
    species : str
        'f_re' for runaway electrons or 'f_hot' for hot electrons
    plot_dir : str
        Directory to save figures
    resonance_params : dict, optional
        If provided, overlay resonance curves. Expected keys:
        - 'omega' : float — wave angular frequency (rad/s)
        - 'k_par' : float — parallel wavenumber (m^-1)
        - 'B0' : float — magnetic field (T)
        - 'n_values' : list of int, optional — harmonics to plot (default [-1, 0, 1])
    """
    
    # Load output data
    print(f"Loading data from {output_file}...")
    do = DREAMOutput(output_file)
    
    # Get distribution function data
    if species == 'f_re':
        f = do.eqsys.f_re.get()
        grid = do.grid.runaway
    elif species == 'f_hot':
        f = do.eqsys.f_hot.get()
        grid = do.grid.hottail
    else:
        raise ValueError("species must be 'f_re' or 'f_hot'")
    
    # Get momentum and pitch grids
    p = grid.p[:]  # momentum in units of m_e*c
    xi = grid.xi[:]  # pitch coordinate (cosine of pitch angle)
    
    # Select radial point (first radial point)
    f_data = f[time_index, 0, :, :]  # [time, radius, xi, p]
    
    # Create meshgrid for plotting
    P, XI = np.meshgrid(p, xi)
    
    # Calculate log(f), handle zeros and negative values
    f_abs = np.abs(f_data)
    f_max = np.max(f_abs)
    if f_max > 0:
        log_f = np.log10(f_abs / f_max + 1e-30)  # Normalize and add small number to avoid log(0)
    else:
        log_f = np.zeros_like(f_abs)
    
    # Get a/R value from grid (matching HDF5 structure: grid/a and grid/R0)
    try:
        a = do.grid.a[:]
        R0 = do.grid.R0[:]
        
        # Handle array or scalar
        a_val = float(a[0]) if hasattr(a, '__len__') and len(a) > 0 else float(a)
        R0_val = float(R0[0]) if hasattr(R0, '__len__') and len(R0) > 0 else float(R0)
        
        a_R = a_val / R0_val
        a_R_text = f'a/R = {a_R:.2f}'
    except Exception as e:
        print(f"Warning: Could not read a/R value: {e}")
        a_R_text = None
    
    # Get toroidal electric field E
    try:
        E_field = do.eqsys.E_field.get()
        E_val = E_field[time_index, 0]
        E_text = f'E = {E_val:.2f} V/m'
    except Exception as e:
        print(f"Warning: Could not read E_field: {e}")
        E_text = None
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Plot contour
    levels = np.linspace(-12, 0, 25)  # Contour levels from -12 to 0
    contour = ax.contourf(P, XI, log_f, levels=levels, cmap='viridis', extend='both')
    
    # Add contour lines
    contour_lines = ax.contour(P, XI, log_f, levels=levels[::2], colors='white', 
                               linewidths=0.5, alpha=0.5)
    
    # Overlay resonance curves if parameters provided
    if resonance_params is not None:
        p_max = np.max(p)
        rl = calculate_resonance_lines(
            p_max,
            omega=resonance_params['omega'],
            k_par=resonance_params['k_par'],
            B0=resonance_params['B0'],
            n_values=resonance_params.get('n_values', None),
        )
        plot_resonance_on_ax(ax, rl, n_values=resonance_params.get('n_values', None))
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)
    
    # Labels and title
    ax.set_xlabel(r'Momentum $p$ [$m_e c$]', fontsize=12)
    ax.set_ylabel(r'Pitch $\xi = \cos(\theta)$', fontsize=12)
    ax.invert_yaxis()  # xi=1 at bottom, xi=-1 at top
    
    time_val = do.grid.t[time_index]
    species_name = 'Runaway electrons ($f_\\mathrm{re}$)' if species == 'f_re' else 'Hot electrons ($f_\\mathrm{hot}$)'
    title = f'{species_name}\nTime = {time_val:.6e} s'
    if a_R_text:
        title += f'\n{a_R_text}'
    if E_text:
        title += f', {E_text}'
    ax.set_title(title, fontsize=14)
    
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(plot_dir, exist_ok=True)
    
    filename = f'logf_contour_{species}_t{time_index}.png'
    save_path = os.path.join(plot_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    
    plt.show()
    
    return fig, ax


def plot_logf_contour_zoomed(output_file='../outputs/dreicer_with_fre_output.h5',
                             time_index=-1, species='f_re', plot_dir='../figures',
                             resonance_params=None,
                             xi_range=(1.0, 0.25), p_range=(0, 30)):
    """
    Plot log(f) contour in p-xi space with custom axis limits (zoomed view).
    
    Parameters:
    -----------
    output_file, time_index, species, plot_dir, resonance_params :
        Same as plot_logf_contour().
    xi_range : tuple
        (xi_max, xi_min) — pitch range for the zoomed view (default (1, 0.25)).
    p_range : tuple
        (p_min, p_max) — momentum range for the zoomed view (default (0, 30)).
    """
    # Load output data
    print(f"Loading data from {output_file}...")
    do = DREAMOutput(output_file)
    
    # Get distribution function data
    if species == 'f_re':
        f = do.eqsys.f_re.get()
        grid = do.grid.runaway
    elif species == 'f_hot':
        f = do.eqsys.f_hot.get()
        grid = do.grid.hottail
    else:
        raise ValueError("species must be 'f_re' or 'f_hot'")
    
    # Get momentum and pitch grids
    p = grid.p[:]  # momentum in units of m_e*c
    xi = grid.xi[:]  # pitch coordinate (cosine of pitch angle)
    
    # Select radial point (first radial point)
    f_data = f[time_index, 0, :, :]  # [time, radius, xi, p]
    
    # Create meshgrid for plotting
    P, XI = np.meshgrid(p, xi)
    
    # Calculate log(f), handle zeros and negative values
    f_abs = np.abs(f_data)
    f_max = np.max(f_abs)
    if f_max > 0:
        log_f = np.log10(f_abs / f_max + 1e-30)
    else:
        log_f = np.zeros_like(f_abs)
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    
    # Plot contour
    levels = np.linspace(-12, 0, 25)
    contour = ax.contourf(P, XI, log_f, levels=levels, cmap='turbo', extend='both')
    
    # Add contour lines
    contour_lines = ax.contour(P, XI, log_f, levels=levels[::2], colors='white',
                               linewidths=0.5, alpha=0.5)
    
    # Overlay resonance curves if parameters provided
    if resonance_params is not None:
        p_max = p_range[1]
        rl = calculate_resonance_lines(
            p_max,
            omega=resonance_params['omega'],
            k_par=resonance_params['k_par'],
            B0=resonance_params['B0'],
            n_values=resonance_params.get('n_values', None),
        )
        plot_resonance_on_ax(ax, rl, n_values=resonance_params.get('n_values', None))
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)
    
    # Labels and title
    ax.set_xlabel(r'Momentum $p$ [$m_e c$]', fontsize=12)
    ax.set_ylabel(r'Pitch $\xi = \cos(\theta)$', fontsize=12)
    
    # Apply zoom limits (ξ=1 at bottom, ξ=0.25 at top)
    ax.set_xlim(p_range[0], p_range[1])
    ax.set_ylim(xi_range[0], xi_range[1])
    
    time_val = do.grid.t[time_index]
    species_name = 'Runaway electrons ($f_\\mathrm{re}$)' if species == 'f_re' else 'Hot electrons ($f_\\mathrm{hot}$)'
    ax.set_title(f'{species_name} (zoomed)\nTime = {time_val:.6e} s', fontsize=14)
    
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(plot_dir, exist_ok=True)
    
    filename = f'logf_contour_{species}_t{time_index}_zoomed.png'
    save_path = os.path.join(plot_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    
    plt.show()
    
    return fig, ax

def plot_multiple_times(output_file='../outputs/dreicer_with_fre_output.h5', 
                       time_indices=None, species='f_re', plot_dir='../figures',
                       resonance_params=None):
    """
    Plot log(f) contour at multiple time steps
    
    Parameters:
    -----------
    output_file : str
        Path to the DREAM output file
    time_indices : list
        List of time indices to plot (None for automatic selection)
    species : str
        'f_re' for runaway electrons or 'f_hot' for hot electrons
    plot_dir : str
        Directory to save figures
    resonance_params : dict, optional
        If provided, overlay resonance curves. Same keys as plot_logf_contour.
    """
    
    # Load output data
    print(f"Loading data from {output_file}...")
    do = DREAMOutput(output_file)
    
    nt = len(do.grid.t[:])
    
    # Select time indices if not provided
    if time_indices is None:
        # Select 5 evenly spaced time points
        time_indices = np.linspace(0, nt-1, 5, dtype=int)
    
    # Get distribution function data
    if species == 'f_re':
        f = do.eqsys.f_re.get()
        grid = do.grid.runaway
    elif species == 'f_hot':
        f = do.eqsys.f_hot.get()
        grid = do.grid.hottail
    else:
        raise ValueError("species must be 'f_re' or 'f_hot'")
    
    # Get momentum and pitch grids
    p = grid.p[:]
    xi = grid.xi[:]
    P, XI = np.meshgrid(p, xi)
    
    # Pre-calculate resonance lines if parameters provided
    resonance_lines = None
    if resonance_params is not None:
        p_max = np.max(p)
        resonance_lines = calculate_resonance_lines(
            p_max,
            omega=resonance_params['omega'],
            k_par=resonance_params['k_par'],
            B0=resonance_params['B0'],
            n_values=resonance_params.get('n_values', None),
        )
    
    # Create subplots
    n_plots = len(time_indices)
    fig, axes = plt.subplots(1, n_plots, figsize=(5*n_plots, 5))
    if n_plots == 1:
        axes = [axes]
    
    for idx, time_idx in enumerate(time_indices):
        # Select radial point (first radial point)
        f_data = f[time_idx, 0, :, :]
        
        # Calculate log(f)
        f_abs = np.abs(f_data)
        f_max = np.max(f_abs)
        if f_max > 0:
            log_f = np.log10(f_abs / f_max + 1e-30)
        else:
            log_f = np.zeros_like(f_abs)
        
        # Plot contour
        levels = np.linspace(-12, 0, 25)
        contour = axes[idx].contourf(P, XI, log_f, levels=levels, cmap='viridis', extend='both')
        
        # Labels
        axes[idx].set_xlabel(r'$p$ [$m_e c$]', fontsize=10)
        if idx == 0:
            axes[idx].set_ylabel(r'$\xi$', fontsize=10)
        axes[idx].invert_yaxis()  # xi=1 at bottom, xi=-1 at top
        
        time_val = do.grid.t[time_idx]
        axes[idx].set_title(f't = {time_val:.6e} s', fontsize=11)
        axes[idx].grid(True, linestyle='--', alpha=0.3)
        
        # Overlay resonance curves on each subplot
        if resonance_lines is not None:
            plot_resonance_on_ax(axes[idx], resonance_lines,
                                 n_values=resonance_params.get('n_values', None))
            axes[idx].legend(loc='upper right', fontsize=8)
    
    # Add colorbar
    cbar = fig.colorbar(contour, ax=axes.tolist(), fraction=0.02, pad=0.04)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)
    
    species_name = 'Runaway electrons' if species == 'f_re' else 'Hot electrons'
    fig.suptitle(f'{species_name} distribution evolution', fontsize=14, y=1.02)
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs(plot_dir, exist_ok=True)
    
    filename = f'logf_evolution_{species}.png'
    save_path = os.path.join(plot_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    
    plt.show()
    
    return fig, axes


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Plot log(f) contour in p-xi space')
    parser.add_argument('--data_dir', type=str, default=None,
                        help='Path to the DREAM output HDF5 file')
    parser.add_argument('--plot_dir', type=str, default='../figures',
                        help='Directory to save plots')

    # Resonance overlay arguments
    parser.add_argument('--omega', type=float, default=2*np.pi*476e6,
                        help='Wave angular frequency (rad/s), default 2π×476 MHz')
    parser.add_argument('--k_par', type=float, default=-41.0,
                        help='Parallel wavenumber (m^-1), default -41')
    parser.add_argument('--B0', type=float, default=1.4,
                        help='Magnetic field (T), default 1.4')
    parser.add_argument('--harmonics', type=str, default='-1,0,1,2',
                        help='Comma-separated harmonic numbers to plot (default: -1,0,1,2)')
    parser.add_argument('--no-resonance', action='store_true',
                        help='Disable resonance curve overlay')
    parser.add_argument('--zoomed', action='store_true',
                        help='Generate zoomed view with custom p/xi limits')
    parser.add_argument('--xi-max', type=float, default=1.0,
                        help='Maximum xi for zoomed view (default: 1.0)')
    parser.add_argument('--xi-min', type=float, default=0.25,
                        help='Minimum xi for zoomed view (default: 0.25)')
    parser.add_argument('--p-max', type=float, default=30,
                        help='Maximum p for zoomed view (default: 30)')
    args = parser.parse_args()

    print("="*60)
    print("DREAM Distribution Function Visualization")
    print("="*60)
    
    # Use provided data_dir or default to whistler output
    output_file = args.data_dir if args.data_dir else '../outputs/quasilinear_whistler_output.h5'
    
    if not Path(output_file).exists():
        print(f"Error: Output file {output_file} not found!")
        print("Please run generate_with_fre.py first.")
        return
    
    # Build resonance parameters dict (unless disabled)
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

    # Plot single time step (last time step)
    print("\n1. Plotting log(f) contour at last time step...")
    plot_logf_contour(output_file, time_index=-1, species='f_re',
                      plot_dir=args.plot_dir, resonance_params=resonance_params)

    # Plot evolution at multiple times
    print("\n2. Plotting log(f) evolution at multiple time steps...")
    plot_multiple_times(output_file, species='f_re',
                        plot_dir=args.plot_dir, resonance_params=resonance_params)
    
    # Zoomed view (if requested)
    if args.zoomed:
        print("\n3. Plotting zoomed log(f) contour...")
        # Only show n=+1 and n=+2 in zoomed view (n=0 and n=-1 not visible)
        zoom_resonance_params = dict(resonance_params) if resonance_params else None
        if zoom_resonance_params is not None:
            zoom_resonance_params['n_values'] = [1, 2]
        plot_logf_contour_zoomed(
            output_file, time_index=-1, species='f_re',
            plot_dir=args.plot_dir, resonance_params=zoom_resonance_params,
            xi_range=(args.xi_max, args.xi_min), p_range=(0, args.p_max)
        )
    
    print("\nVisualization completed!")


if __name__ == "__main__":
    main()
