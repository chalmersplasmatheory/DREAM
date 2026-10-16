#!/usr/bin/env python3
#
# Plots the current densities, the temperatures and densities for electrons and ions.
#
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

from DREAM.DREAMOutput import DREAMOutput

LOADFILE = "output.h5"

def plotCurrents(do, ax=None):
    """
    Plots the profiles for the Ohmic, Bootstrap and total currents, given a DREAMOutput object.
    """
    if ax is None:
        ax = plt.axes()

    ax.plot(do.eqsys.j_ohm.radius, do.eqsys.j_ohm.data[-1,:]*1e-6, label="Ohmic")
    ax.plot(do.eqsys.j_bs.radius,  do.eqsys.j_bs.data[-1,:]*1e-6,  label="Bootstrap")
    ax.plot(do.eqsys.j_tot.radius, do.eqsys.j_tot.data[-1,:]*1e-6, label="Total")

    ax.set_xlabel("minor radius [m]")
    ax.set_ylabel(r"current density [MA m$^{-2}$]")
    ax.legend()

    return ax

def plotTemperatures(do, ax=None):
    """
    Plots the temperature profiles for the electron and ions, given a DREAMOutput object.
    """
    if ax is None:
        ax = plt.axes()

    ax.plot(do.eqsys.T_cold.radius, do.eqsys.T_cold.data[-1,:]*1e-3, label="Electrons")
    for N_i, W_i, name in zip(do.eqsys.N_i.data, do.eqsys.W_i.data, do.eqsys.W_i.ions.names):
        T_i = W_i / (1.5 * N_i * sc.e)
        ax.plot(do.eqsys.W_i.grid.r, T_i[-1,:]*1e-3, label=name)

    ax.set_xlabel("minor radius [m]")
    ax.set_ylabel("temperature [keV]")
    ax.legend()

    return ax

def plotDensities(do, ax=None):
    """
    Plots the density profiles for the electron and ions, given a DREAMOutput object.
    """
    if ax is None:
        ax = plt.axes()

    ax.plot(do.eqsys.n_cold.radius, do.eqsys.n_cold.data[-1,:]*1e-20, label="Electrons")
    for N_i, name in zip(do.eqsys.N_i.data, do.eqsys.N_i.ions.names):
        ax.plot(do.eqsys.N_i.grid.r, N_i[-1,:]*1e-20, label=name)

    ax.set_xlabel("minor radius [m]")
    ax.set_ylabel(r"density [$10^{20}$ m$^{-3}$]")
    ax.legend()

    return ax

if __name__ == '__main__':

    do = DREAMOutput(LOADFILE)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 10), sharex=True)
    plotCurrents(do, ax1)
    plotTemperatures(do, ax2)
    plotDensities(do, ax3)

    plt.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.show()
