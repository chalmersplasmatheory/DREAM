#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')
from DREAM import DREAMOutput

def plot_distribution_contours(do, time_index, title_suffix, filename):
    """
    Plot contours of the distribution function at a specific time
    """
    # Get the distribution function at the specified time
    f_hot = do.eqsys.f_hot[time_index, 0, :, :]  # t, r, xi, p
    
    # Normalize by the maximum value
    f_max = np.max(f_hot)
    if f_max > 0:
        F = f_hot / f_max
    else:
        F = f_hot
    
    # Take logarithm (avoid log(0) by adding a small value)
    F_log = np.log10(F + 1e-10)
    
    # Get momentum and pitch grids
    p = do.grid.hottail.p
    xi = do.grid.hottail.xi
    
    # Create meshgrid for plotting
    P, XI = np.meshgrid(p, xi, indexing='ij')
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    cont = ax.contourf(P, XI, F_log.T, levels=50, cmap='plasma')
    ax.set_xlabel(r'Momentum $p/m_ec$')
    ax.set_ylabel(r'Pitch $\xi$')
    ax.set_title(f'Log10 of distribution function F = f/max(f) at {title_suffix}')
    plt.colorbar(cont, ax=ax, label=r'$\log_{10}(F)$')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def main():
    # Load the output
    do = DREAMOutput('../outputs/spitzer_output.h5')
    
    # Get time grid and tauEETh
    t = do.grid.t[:]
    tauEETh_val = 5.204e-5  # From our calculation
    
    # Calculate thermal collision times elapsed at each time point
    thermal_times = t / tauEETh_val
    
    print(f"Simulation time points: {len(t)}")
    print(f"Final time in thermal collision times: {thermal_times[-1]:.1f}")
    
    # Find time indices corresponding to approximately 0, 167, 333, and 500 thermal collision times
    target_times = [0, 167, 333, 500]
    time_indices = []
    
    for target in target_times:
        # Find the closest time index
        idx = np.argmin(np.abs(thermal_times - target))
        time_indices.append(idx)
        print(f"Target tau_th = {target}: actual = {thermal_times[idx]:.1f} at time index {idx}")
    
    # Plot distributions at these times (with log scale)
    for i, (target_time, time_idx) in enumerate(zip(target_times, time_indices)):
        plot_distribution_contours(
            do, 
            time_idx, 
            f"τ_th = {target_time}", 
            f"../figures/spitzer_distribution_log_t{target_time}.png"
        )
    
    # Also plot the time evolution of some key quantities
    plt.figure(figsize=(10, 6))
    j_tot = do.eqsys.j_tot[:, 0]
    plt.plot(thermal_times, j_tot, 'b-', linewidth=2)
    plt.xlabel(r'Thermal collision times $\tau_{th}$')
    plt.ylabel('Total current density (A/m²)')
    plt.title('Total Current Density Evolution')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('../figures/spitzer_current_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Plots saved to ../figures/")

if __name__ == "__main__":
    main()