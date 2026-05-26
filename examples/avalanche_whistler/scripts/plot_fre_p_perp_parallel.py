#!/usr/bin/env python3
#
# Plot p_perp vs p_parallel for the runaway electron distribution (f_re)
#

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput

import argparse
import os

parser = argparse.ArgumentParser(description='Plot p_perp vs p_parallel for runaway electrons')
parser.add_argument('--data_dir', type=str, default='../outputs/dreicer_with_fre_output.h5',
                    help='Path to the DREAM output HDF5 file')
parser.add_argument('--plot_dir', type=str, default='../figures',
                    help='Directory to save plots')
args = parser.parse_args()

os.makedirs(args.plot_dir, exist_ok=True)

# Load data
do = DREAMOutput(args.data_dir)

print("Available quantities in eqsys:")
print([attr for attr in dir(do.eqsys) if not attr.startswith('_')])

print("\nChecking for f_re:")
if hasattr(do.eqsys, 'f_re'):
    print("f_re exists")
    f_re_shape = do.eqsys.f_re[:].shape
    print(f"f_re shape: {f_re_shape}")
else:
    print("f_re does not exist")
    do.close()
    exit()

# Get momentum grid for runaway electrons
p_re = do.grid.runaway.p[:]
xi_re = do.grid.runaway.xi[:]
print(f"\nRunaway grid: {len(p_re)} p-points, {len(xi_re)} xi-points")

# Calculate p_parallel and p_perp for runaway electrons
# For each (p, xi) pair, we have:
# p_parallel = p * xi
# p_perp = p * sqrt(1 - xi^2)
p_mesh, xi_mesh = np.meshgrid(p_re, xi_re, indexing='ij')  # Shape: (np, nxi)
p_parallel = p_mesh * xi_mesh
p_perp = p_mesh * np.sqrt(1 - xi_mesh**2)

# Get distribution function at last time step
time_index = -1
f_re = do.eqsys.f_re[time_index, 0, :, :]  # Shape: (nxi, np)

# Transpose to match the dimensions of p_parallel and p_perp
f_re_transposed = f_re.T  # Now shape is (np, nxi)

# Avoid taking log of zero or negative values
f_re_positive = np.where(f_re_transposed > 1e-30, f_re_transposed, 1e-30)

# Create figure with scatter plot
plt.figure(figsize=(10, 8))

# Flatten arrays for scatter plot
p_parallel_flat = p_parallel.flatten()
p_perp_flat = p_perp.flatten()
f_re_flat = f_re_positive.flatten()

# Filter out points with very low distribution values for clarity
mask = f_re_flat > 1e-20
p_parallel_filtered = p_parallel_flat[mask]
p_perp_filtered = p_perp_flat[mask]
f_re_filtered = f_re_flat[mask]

# Create scatter plot with color representing log10 of distribution function
scatter = plt.scatter(p_parallel_filtered, p_perp_filtered, 
                      c=np.log10(f_re_filtered), 
                      cmap='plasma', s=1, alpha=0.8)
cbar = plt.colorbar(scatter, label=r'$\log_{10}(f_{\rm re})$')
cbar.ax.tick_params(labelsize=12)

plt.xlabel(r'$p_\parallel$ ($m_e c$)', fontsize=14)
plt.ylabel(r'$p_\perp$ ($m_e c$)', fontsize=14)
plt.title('Runaway electron distribution in momentum space\n(last time step)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.7)

# Set axis limits for better visualization
plt.xlim(-10.0, 10.0)
plt.ylim(0, 10.0)

# Save figure
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_fre_p_perp_parallel_scatter.png'), dpi=300)
plt.show()

# Create contour plot version
plt.figure(figsize=(10, 8))

# Determine appropriate contour levels
min_val = np.log10(np.min(f_re_positive))
max_val = np.log10(np.max(f_re_positive))
print(f"Min log10(f_re): {min_val:.2f}")
print(f"Max log10(f_re): {max_val:.2f}")

# Ensure we have valid contour levels
if max_val > min_val:
    levels = np.linspace(min_val, max_val, 30)
    contour = plt.contourf(p_parallel, p_perp, np.log10(f_re_positive), 
                           levels=levels, cmap='plasma', extend='both')
    cbar = plt.colorbar(contour, label=r'$\log_{10}(f_{\rm re})$')
    cbar.ax.tick_params(labelsize=12)
else:
    # If all values are the same, just plot the data
    contour = plt.contourf(p_parallel, p_perp, f_re_positive, cmap='plasma')
    cbar = plt.colorbar(contour, label=r'$f_{\rm re}$')
    cbar.ax.tick_params(labelsize=12)

plt.xlabel(r'$p_\parallel$ ($m_e c$)', fontsize=14)
plt.ylabel(r'$p_\perp$ ($m_e c$)', fontsize=14)
plt.title('Runaway electron distribution in momentum space\n(last time step, contour)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.7)

# Set axis limits
plt.xlim(-10.0, 10.0)
plt.ylim(0, 10.0)

# Save figure
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_fre_p_perp_parallel_contour.png'), dpi=300)
plt.show()

# Create another version showing multiple time steps
plt.figure(figsize=(15, 5))

# Select three time steps to compare
time_steps = [0, len(do.grid.t)//2, -1]  # Initial, middle, final
time_labels = ['Initial', 'Middle', 'Final']

for i, (t_idx, label) in enumerate(zip(time_steps, time_labels)):
    plt.subplot(1, 3, i+1)
    
    # Get distribution function at this time step
    f_re_t = do.eqsys.f_re[t_idx, 0, :, :].T  # Transpose to (np, nxi)
    f_re_t_positive = np.where(f_re_t > 1e-30, f_re_t, 1e-30)
    
    # Determine appropriate contour levels for this time step
    min_val = np.log10(np.min(f_re_t_positive))
    max_val = np.log10(np.max(f_re_t_positive))
    
    plt.xlabel(r'$p_\parallel$ ($m_e c$)', fontsize=12)
    plt.ylabel(r'$p_\perp$ ($m_e c$)', fontsize=12)
    plt.title(f'{label} (t = {do.grid.t[t_idx]*1e6:.1f} μs)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(-10.0, 10.0)
    plt.ylim(0, 10.0)
    
    # Only plot contours if we have variation in values
    if max_val > min_val:
        levels = np.linspace(min_val, max_val, 20)
        contour = plt.contourf(p_parallel, p_perp, np.log10(f_re_t_positive), 
                               levels=levels, cmap='plasma', extend='both')
    else:
        plt.contourf(p_parallel, p_perp, f_re_t_positive, cmap='plasma')

# Add a colorbar for all subplots
plt.tight_layout()
plt.savefig(os.path.join(args.plot_dir, 'dreicer_fre_p_perp_parallel_evolution.png'), dpi=300)
plt.show()

# Print some statistics
print("Momentum space visualization:")
print(f"Runaway grid: {len(p_re)} p-points, {len(xi_re)} xi-points")
print(f"Max p_parallel: {np.max(p_parallel):.2f} m_e c")
print(f"Max p_perp: {np.max(p_perp):.2f} m_e c")

# Check runaway electron density
n_re = do.eqsys.n_re[:, 0]
print(f"\nRunaway electron density:")
print(f"  Initial: {n_re[0]:.2e} m^-3")
print(f"  Final:   {n_re[-1]:.2e} m^-3")

# Close the file
do.close()