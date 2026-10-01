#!/usr/bin/env python3
#
# Plot p_perp vs p_parallel for the runaway electron distribution (f_re) 
# from the Spitzer case
#

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import DREAMOutput

def plot_fre_p_perp_parallel(do, time_index, title_suffix, filename):
    """
    Plot contours of the runaway electron distribution function in p_perp-p_parallel coordinates
    with p_parallel as x-axis and p_perp as y-axis
    """
    # Check if runaway grid exists
    if not hasattr(do.grid, 'runaway') or do.grid.runaway is None:
        print("No runaway electron grid found in the data")
        return False
    
    # Get the runaway distribution function at the specified time
    f_re = do.eqsys.f_re[time_index, 0, :, :]  # t, r, xi, p
    
    # Transpose to match the dimensions of p_parallel and p_perp
    f_re_transposed = f_re.T  # Now shape is (np, nxi)
    
    # Avoid taking log of zero or negative values
    f_re_positive = np.where(f_re_transposed > 1e-30, f_re_transposed, 1e-30)
    
    # Take logarithm
    F_re_log = np.log10(f_re_positive)
    
    # Get momentum and pitch grids
    p_re = do.grid.runaway.p[:]
    xi_re = do.grid.runaway.xi[:]
    
    # Create 2D meshes for p and xi
    P, XI = np.meshgrid(p_re, xi_re, indexing='ij')
    
    # Convert to p_perp and p_parallel
    P_PERP = P * np.sqrt(1 - XI**2)
    P_PARALLEL = P * XI
    
    # Plot with p_parallel as x-axis and p_perp as y-axis
    fig, ax = plt.subplots(figsize=(10, 8))
    cont = ax.contourf(P_PARALLEL, P_PERP, F_re_log, levels=50, cmap='plasma')
    ax.set_xlabel(r'$p_\parallel/m_ec$')
    ax.set_ylabel(r'$p_\perp/m_ec$')
    ax.set_title(f'Log10 of runaway electron distribution function f_re at {title_suffix}')
    plt.colorbar(cont, ax=ax, label=r'$\log_{10}(f_{re})$')
    
    # Set equal aspect ratio for better visualization
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    return True

def main():
    # Load the output
    try:
        do = DREAMOutput('../outputs/spitzer_output.h5')
    except Exception as e:
        print(f"Error loading output file: {e}")
        return
    
    # Check if runaway grid exists
    if not hasattr(do.grid, 'runaway') or do.grid.runaway is None:
        print("No runaway electron grid found in the data.")
        do.close()
        return
    
    # Get time grid
    t = do.grid.t[:]
    
    print(f"Simulation time points: {len(t)}")
    print(f"Runaway grid: {len(do.grid.runaway.p)} p-points, {len(do.grid.runaway.xi)} xi-points")
    
    # Plot distribution at final time
    success = plot_fre_p_perp_parallel(
        do, 
        -1, 
        f"final time", 
        f"../figures/spitzer_fre_p_parallel_perp_log_t{len(t)-1}.png"
    )
    
    if success:
        print("Plot saved to ../figures/")
    else:
        print("Failed to generate plot")
    
    do.close()

if __name__ == "__main__":
    main()