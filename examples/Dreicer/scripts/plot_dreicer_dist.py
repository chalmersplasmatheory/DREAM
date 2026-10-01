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


do = DREAMOutput('../outputs/dreicer_with_fre_output.h5')

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
    plt.savefig('../figures/dreicer_f_hot_2d.png', dpi=300)
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
plt.savefig('../figures/dreicer_density_evolution.png', dpi=300)
plt.show()

# Plot runaway electron generation rate
plt.figure(figsize=(10, 6))

# Calculate numerical derivative of n_re
dt = np.diff(t)
dn_re = np.diff(n_re)
gamma_re = np.zeros_like(t)
gamma_re[1:] = dn_re/(dt*n_re[1:])  # Growth rate

# Mask invalid values
mask = (n_re > 1e5) & (~np.isnan(gamma_re)) & (~np.isinf(gamma_re))
plt.plot(t[mask], gamma_re[mask], 'r-', linewidth=2)
plt.xlabel(r'Time (s)')
plt.ylabel(r'Dreicer growth rate $\gamma_{\rm Dreicer}$ (s$^{-1}$)')
plt.title('Dreicer Growth Rate Evolution')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('../figures/dreicer_growth_rate.png', dpi=300)
plt.show()

# Print some statistics
print("Simulation Statistics:")
print(f"Final runaway electron density: {n_re[-1]:.2e} m^-3")
print(f"Total electrons: {n_cold[0] + n_hot[0] + n_re[0]:.2e} m^-3 (initial)")
print(f"Total electrons: {n_cold[-1] + n_hot[-1] + n_re[-1]:.2e} m^-3 (final)")

# Close the DREAMOutput file
do.close()