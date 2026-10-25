#!/usr/bin/env python3
#
# Plot p_perp vs p_parallel for the hot electron distribution
#

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput

import argparse
import os

parser = argparse.ArgumentParser(description='Plot p_perp vs p_parallel for hot electrons')
parser.add_argument('--data_dir', type=str, default='../outputs/dreicer_with_fre_output.h5',
                    help='Path to the DREAM output HDF5 file')
parser.add_argument('--plot_dir', type=str, default='../figures',
                    help='Directory to save plots')
args = parser.parse_args()

os.makedirs(args.plot_dir, exist_ok=True)

# Load data
do = DREAMOutput(args.data_dir)

# Get the last time step
time_index = -1

# Get momentum grid
p = do.grid.hottail.p[:]
xi = do.grid.hottail.xi[:]

# Calculate p_parallel and p_perp
# For each (p, xi) pair, we have:
# p_parallel = p * xi
# p_perp = p * sqrt(1 - xi^2)
p_mesh, xi_mesh = np.meshgrid(p, xi, indexing='ij')  # Shape: (np, nxi)
p_parallel = p_mesh * xi_mesh
p_perp = p_mesh * np.sqrt(1 - xi_mesh**2)

# Get distribution function at last time step
f_hot = do.eqsys.f_hot[time_index, 0, :, :]  # Shape: (nxi, np)

# Transpose to match the dimensions of p_parallel and p_perp
f_hot_transposed = f_hot.T  # Now shape is (np, nxi)

# Avoid taking log of zero or negative values
f_hot_positive = np.where(f_hot_transposed > 1e-30, f_hot_transposed, 1e-30)

# Create figure with scatter plot
plt.figure(figsize=(10, 8))

# Flatten arrays for scatter plot
p_parallel_flat = p_parallel.flatten()
p_perp_flat = p_perp.flatten()
f_hot_flat = f_hot_positive.flatten()

# Create scatter plot with color representing log10 of distribution function
scatter = plt.scatter(p_parallel_flat, p_perp_flat, 
                      c=np.log10(f_hot_flat), 
                      cmap='plasma', s=1, alpha=0.8)
cbar = plt.colorbar(scatter, label=r'$\log_{10}(f_{\rm hot})$')
cbar.ax.tick_params(labelsize=12)

plt.xlabel(r'$p_\parallel$ ($m_e c$)', fontsize=14)
plt.ylabel(r'$p_\perp$ ($m_e c$)', fontsize=14)
plt.title('Hot electron distribution in momentum space\n(last time step)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.7)

# Set axis limits for better visualization
plt.xlim(-1.0, 1.0)
plt.ylim(0, 1.1)

# Save figure
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_p_perp_parallel_scatter.png'), dpi=300)
plt.show()

# Create contour plot version
plt.figure(figsize=(10, 8))

# Use p_parallel and p_perp directly for plotting
levels = np.linspace(np.log10(np.min(f_hot_positive)), 
                     np.log10(np.max(f_hot_positive)), 30)

contour = plt.contourf(p_parallel, p_perp, np.log10(f_hot_positive), 
                       levels=levels, cmap='plasma', extend='both')
cbar = plt.colorbar(contour, label=r'$\log_{10}(f_{\rm hot})$')
cbar.ax.tick_params(labelsize=12)

plt.xlabel(r'$p_\parallel$ ($m_e c$)', fontsize=14)
plt.ylabel(r'$p_\perp$ ($m_e c$)', fontsize=14)
plt.title('Hot electron distribution in momentum space\n(last time step, contour)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.7)

# Set axis limits
plt.xlim(-1.0, 1.0)
plt.ylim(0, 1.1)

# Save figure
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_p_perp_parallel_contour.png'), dpi=300)
plt.show()

# Create another version showing multiple time steps
plt.figure(figsize=(15, 5))

# Select three time steps to compare
time_steps = [0, len(do.grid.t)//2, -1]  # Initial, middle, final
time_labels = ['Initial', 'Middle', 'Final']

for i, (t_idx, label) in enumerate(zip(time_steps, time_labels)):
    plt.subplot(1, 3, i+1)
    
    # Get distribution function at this time step
    f_hot_t = do.eqsys.f_hot[t_idx, 0, :, :].T  # Transpose to (np, nxi)
    f_hot_t_positive = np.where(f_hot_t > 1e-30, f_hot_t, 1e-30)
    
    # Create contour plot
    levels = np.linspace(np.log10(np.min(f_hot_t_positive)), 
                         np.log10(np.max(f_hot_t_positive)), 20)
    
    contour = plt.contourf(p_parallel, p_perp, np.log10(f_hot_t_positive), 
                           levels=levels, cmap='plasma', extend='both')
    
    plt.xlabel(r'$p_\parallel$ ($m_e c$)', fontsize=12)
    plt.ylabel(r'$p_\perp$ ($m_e c$)', fontsize=12)
    plt.title(f'{label} (t = {do.grid.t[t_idx]*1e6:.1f} μs)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(-1.0, 1.0)
    plt.ylim(0, 1.1)

# Add a colorbar for all subplots
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_p_perp_parallel_evolution.png'), dpi=300)
plt.show()

# Print some statistics
print("Momentum space visualization:")
print(f"Momentum grid: {len(p)} p-points, {len(xi)} xi-points")
print(f"Max p_parallel: {np.max(p_parallel):.2f} m_e c")
print(f"Max p_perp: {np.max(p_perp):.2f} m_e c")

# Close the file
do.close()