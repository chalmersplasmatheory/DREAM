


#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')
from DREAM import DREAMOutput

def plot_p_perp_parallel(do, time_index, title_suffix, filename):
    """
    Plot contours of the distribution function in p_perp-p_parallel coordinates
    with p_parallel as x-axis and p_perp as y-axis
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
    
    # Truncate values above -1
    F_log_truncated = np.where(F_log > -1, -1, F_log)
    
    # Get momentum and pitch grids
    p = do.grid.hottail.p[:]
    xi = do.grid.hottail.xi[:]
    
    # Create 2D meshes for p and xi
    P, XI = np.meshgrid(p, xi, indexing='ij')
    
    # Convert to p_perp and p_parallel
    P_PERP = P * np.sqrt(1 - XI**2)
    P_PARALLEL = P * XI
    
    # Plot with p_parallel as x-axis and p_perp as y-axis
    fig, ax = plt.subplots(figsize=(10, 8))
    cont = ax.contourf(P_PARALLEL, P_PERP, F_log_truncated.T, levels=20, cmap='plasma', vmin=-10, vmax=-1)
    ax.set_xlabel(r'$p_\parallel/m_ec$')
    ax.set_ylabel(r'$p_\perp/m_ec$')
    ax.set_title(f'Log10 of distribution function F = f/max(f) at {title_suffix}')
    plt.colorbar(cont, ax=ax, label=r'$\log_{10}(F)$', extend='max')
    
    # Set equal aspect ratio for better visualization
    ax.set_aspect('equal', adjustable='box')
    
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
    
    # Plot distributions at these times in p_parallel-p_perp coordinates (with log scale)
    for i, (target_time, time_idx) in enumerate(zip(target_times, time_indices)):
        plot_p_perp_parallel(
            do, 
            time_idx, 
            f"τ_th = {target_time}", 
            f"../figures/spitzer_distribution_p_parallel_perp_log_truncated_t{target_time}.png"
        )
    
    print("Plots saved to ../figures/")

if __name__ == "__main__":
    main()