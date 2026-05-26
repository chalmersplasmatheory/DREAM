#!/usr/bin/env python3
#
# Plot Dreicer generation results
# ######

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput


import argparse
import os

parser = argparse.ArgumentParser(description='Plot Dreicer generation results')
parser.add_argument('--data_dir', type=str, default='../outputs/dreicer_with_fre_output.h5',
                    help='Path to the DREAM output HDF5 file')
parser.add_argument('--plot_dir', type=str, default='../figures',
                    help='Directory to save plots')
args = parser.parse_args()

os.makedirs(args.plot_dir, exist_ok=True)

do = DREAMOutput(args.data_dir)

timeindex = -1
# For Dreicer case, we only have f_hot since runaway grid is disabled
if hasattr(do.eqsys, 'f_hot'):
    fhot = do.eqsys.f_hot[timeindex, 0, :, :]  # Get last time slice, first radial grid point
    
    # Plot 2D distribution function
    plt.figure(figsize=(10, 6))
    # Transpose correctly to match dimensions
    fhot_plot = fhot
    # Avoid taking log of zero or negative values
    fhot_plot = np.where(fhot_plot > 0, fhot_plot, 1e-30)
    
    # Create proper meshgrid for contour plot
    p = do.grid.hottail.p[:]
    xi = do.grid.hottail.xi[:]
    P, XI = np.meshgrid(p, xi)
    
    levels = np.linspace(-20, np.log10(np.max(fhot_plot)), 50)
    plt.contourf(P, XI, np.log10(fhot_plot), levels=levels, cmap='plasma')
    plt.colorbar(label=r'$\log_{10}(f_{\rm hot})$')
    plt.xlabel(r'Momentum $p/m_ec$')
    plt.ylabel(r'Pitch $\xi$')
    plt.title(r'Hot electron distribution function (log scale)')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'dreicer_f_hot_2d.png'), dpi=300)
    plt.show()

# Plot electron density evolution
plt.figure(figsize=(10, 6))

t = do.grid.t[:]
n_cold = do.eqsys.n_cold[:, 0]
n_hot = do.eqsys.n_hot[:, 0] 
n_re = do.eqsys.n_re[:, 0]

plt.plot(t, n_cold, label=r'$n_{\rm cold}$')
plt.plot(t, n_hot, label=r'$n_{\rm hot}$')
plt.plot(t, n_re, label=r'$n_{\rm re}$')
plt.plot(t, n_cold + n_hot + n_re, '--', label=r'$n_{\rm total}$')

plt.xlabel(r'Time (s)')
plt.ylabel(r'Electron density (m$^{-3}$)')
plt.title('Electron density evolution')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.yscale('log')
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_density_evolution.png'), dpi=300)
plt.show()

# Plot runaway electron generation rate
plt.figure(figsize=(10, 6))

# Calculate numerical derivative of n_re
dt = np.diff(t)
dn_re = np.diff(n_re)
gamma_re = np.zeros_like(t)
gamma_re[1:] = dn_re/(dt*n_re[1:])  # Growth rate

# Mask invalid values
t_avg = 0.5 * (t[:-1] + t[1:])  # Average time for growth rate
gamma_re_valid = gamma_re[1:]
mask = (n_re[1:] > 1e5) & (~np.isnan(gamma_re_valid)) & (~np.isinf(gamma_re_valid))

# Skip first 0.1 seconds (transient phase)
time_cutoff = 0.1  # seconds
time_mask = t_avg >= time_cutoff
combined_mask = mask & time_mask

# Plot only data after 0.1s
plt.plot(t_avg[combined_mask], gamma_re_valid[combined_mask], 'r-', linewidth=2, label='Growth rate')

# Get the final growth rate
if np.any(combined_mask):
    final_gamma = gamma_re_valid[combined_mask][-1]
    final_time = t_avg[combined_mask][-1]
    
    # Add annotation showing final growth rate
    plt.annotate(f'Final $\gamma$ = {final_gamma:.3e} s$^{{-1}}$\nat t = {final_time:.3f} s',
                xy=(final_time, final_gamma),
                xytext=(0.7*final_time, 0.7*final_gamma),
                fontsize=11,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.3', lw=2))
    
    print(f"\nFinal Growth Rate Statistics:")
    print(f"Time: {final_time:.6f} s")
    print(f"Growth rate: {final_gamma:.6e} s^-1")
    
    # Also calculate average growth rate in the last 10% of simulation
    valid_times = t_avg[combined_mask]
    valid_gammas = gamma_re_valid[combined_mask]
    if len(valid_gammas) > 10:
        last_10_percent_idx = int(0.9 * len(valid_gammas))
        avg_gamma_last = np.mean(valid_gammas[last_10_percent_idx:])
        print(f"Average growth rate (last 10%): {avg_gamma_last:.6e} s^-1")
else:
    print("\nWarning: No valid growth rate data after 0.1s")

plt.xlabel(r'Time (s)')
plt.ylabel(r'Growth rate $\gamma = (1/n_{re}) dn_{re}/dt$ (s$^{-1}$)')
plt.title(f'Total Runaway Electron Growth Rate (t ≥ {time_cutoff} s)')
plt.legend(loc='best')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_growth_rate.png'), dpi=300)
plt.show()

# Print some statistics
print("Simulation Statistics:")
print(f"Final runaway electron density: {n_re[-1]:.2e} m^-3")
print(f"Total electrons: {n_cold[0] + n_hot[0] + n_re[0]:.2e} m^-3 (initial)")
print(f"Total electrons: {n_cold[-1] + n_hot[-1] + n_re[-1]:.2e} m^-3 (final)")

# Close the DREAMOutput file
do.close()