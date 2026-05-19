#!/usr/bin/env python3
#
# Plot log(f) contour in p-xi space for runaway electrons (f_re) and hot electrons (f_hot)
#

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import DREAMOutput

def plot_logf_contour(output_file='../outputs/dreicer_with_fre_output.h5', 
                      time_index=-1, species='f_re'):
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
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Plot contour
    levels = np.linspace(-12, 0, 25)  # Contour levels from -12 to 0
    contour = ax.contourf(P, XI, log_f, levels=levels, cmap='viridis', extend='both')
    
    # Add contour lines
    contour_lines = ax.contour(P, XI, log_f, levels=levels[::2], colors='white', 
                               linewidths=0.5, alpha=0.5)
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)
    
    # Labels and title
    ax.set_xlabel(r'Momentum $p$ [$m_e c$]', fontsize=12)
    ax.set_ylabel(r'Pitch $\xi = \cos(\theta)$', fontsize=12)
    
    time_val = do.grid.t[time_index]
    species_name = 'Runaway electrons ($f_\\mathrm{re}$)' if species == 'f_re' else 'Hot electrons ($f_\\mathrm{hot}$)'
    ax.set_title(f'{species_name}\nTime = {time_val:.6e} s', fontsize=14)
    
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_dir = Path('../figures')
    output_dir.mkdir(exist_ok=True)
    
    filename = f'logf_contour_{species}_t{time_index}.png'
    save_path = output_dir / filename
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    
    plt.show()
    
    return fig, ax


def plot_multiple_times(output_file='../outputs/dreicer_with_fre_output.h5', 
                       time_indices=None, species='f_re'):
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
        
        time_val = do.grid.t[time_idx]
        axes[idx].set_title(f't = {time_val:.6e} s', fontsize=11)
        axes[idx].grid(True, linestyle='--', alpha=0.3)
    
    # Add colorbar
    cbar = fig.colorbar(contour, ax=axes.tolist(), fraction=0.02, pad=0.04)
    cbar.set_label(r'$\log_{10}(f/f_{\mathrm{max}})$', fontsize=12)
    
    species_name = 'Runaway electrons' if species == 'f_re' else 'Hot electrons'
    fig.suptitle(f'{species_name} distribution evolution', fontsize=14, y=1.02)
    
    plt.tight_layout()
    
    # Save figure
    output_dir = Path('../figures')
    output_dir.mkdir(exist_ok=True)
    
    filename = f'logf_evolution_{species}.png'
    save_path = output_dir / filename
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {save_path}")
    
    plt.show()
    
    return fig, axes


def main():
    """Main function"""
    print("="*60)
    print("DREAM Distribution Function Visualization")
    print("="*60)
    
    output_file = '../outputs/dreicer_with_fre_output.h5'
    
    if not Path(output_file).exists():
        print(f"Error: Output file {output_file} not found!")
        print("Please run generate_with_fre.py first.")
        return
    
    # Plot single time step (last time step)
    print("\n1. Plotting log(f) contour at last time step...")
    plot_logf_contour(output_file, time_index=-1, species='f_re')
    
    # Plot evolution at multiple times
    print("\n2. Plotting log(f) evolution at multiple time steps...")
    plot_multiple_times(output_file, species='f_re')
    
    print("\nVisualization completed!")


if __name__ == "__main__":
    main()
