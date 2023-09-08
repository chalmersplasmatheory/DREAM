#!/usr/bin/env python3
#
# This example evaluates the mean RE speed at a fixed electric field and for varying
# densities of neon and deuterium. Similar plots can be found in Benjamin Buchholz'
# MSc thesis (to be released...).
#
# Run as
#
#   $ ./contourplot.py
#
############################################################################

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

import DREAM
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

# physical parameters
TEMPERATURE = 1e2
ELECTRIC_FIELD = 10
MAGNETIC_FIELD = 5.25


def generate_settings(nD, nNe):
    """
    Generate DREAM settings object.
    """
    nr = 1
    if isinstance(nD, (np.ndarray, list)):
        assert isinstance(nNe, float)
        nr = len(nD)
    elif isinstance(nNe, (np.ndarray, list)):
        assert isinstance(nD, float)
        nr = len(nNe)

    ds = DREAMSettings()

    ds.collisions.collfreq_mode         = Collisions.COLLFREQ_MODE_FULL
    ds.collisions.collfreq_type         = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda              = Collisions.LNLAMBDA_ENERGY_DEPENDENT

    ds.radialgrid.setB0(MAGNETIC_FIELD)
    ds.radialgrid.setMinorRadius(1)
    ds.radialgrid.setWallRadius(1)
    ds.radialgrid.setNr(nr)

    ds.timestep.setTmax(1)
    ds.timestep.setNt(1)

    ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ds.eqsys.n_re.setFluidSpeed(Runaways.FLUID_SPEED_MODE_ANALYTICAL_DIST)
    ds.eqsys.n_re.setEceff(Eceff=Runaways.COLLQTY_ECEFF_MODE_FULL)
    ds.eqsys.n_re.setInitialProfile(1e10)

    ds.eqsys.E_field.setPrescribedData(ELECTRIC_FIELD)
    ds.eqsys.T_cold.setPrescribedData(TEMPERATURE)

    t = np.linspace(0, ds.timestep.tmax, ds.timestep.nt)
    dr = ds.radialgrid.a / ( ds.radialgrid.nr + 1 )
    r = np.linspace(dr, ds.radialgrid.nr * dr, ds.radialgrid.nr)

    def density(n0, nz):
        n = np.zeros(shape=(nz+1, ds.timestep.nt, ds.radialgrid.nr))
        n[1,:,:] += n0
        return n

    ds.eqsys.n_i.addIon(name='D', Z=1, n=density(nD, 1), r=r, t=t, iontype=IonSpecies.IONS_PRESCRIBED)
    ds.eqsys.n_i.addIon(name='Ne', Z=10, n=density(nNe, 10), r=r, t=t, iontype=IonSpecies.IONS_PRESCRIBED)

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    ds.solver.setType(Solver.LINEAR_IMPLICIT)

    return ds


def getSpeed(do):
    """
    Returns the mean RE speed u_re/c, obtained by the current density of the form j_re = e*u_re*n_re,
    from given DREAMOutput object.
    """
    jre = do.eqsys.j_re.data[-1,:]
    nre = do.eqsys.n_re.data[-1,:]
    # print(jre)
    return  jre / (nre * sc.e * sc.c)


if __name__ == '__main__':

    nscans = 50
    nD_arr = np.logspace(17.6, 22, nscans)
    nNe_arr = np.logspace(17.6, 22, nscans)

    speeds = np.empty((nscans, nscans))
    for ns, nD in enumerate(nD_arr):
        ds = generate_settings(nD, nNe_arr)
        do = DREAM.runiface(ds, quiet=True)
        speeds[:,ns] = getSpeed(do)
        # print(speeds)

    ax = plt.axes()
    cp = ax.contourf(nD_arr, nNe_arr, speeds, cmap='jet', levels=np.linspace(.82, 1, 50))
    cb = plt.colorbar(mappable=cp, ax=ax)
    ax.set_yscale('log')
    ax.set_xscale('log')

    cb.ax.set_ylabel(r"$u_{\rm re} / c$")
    ax.set_xlabel(r"$\log_{10}\left(n_{D(+1)}\,[{\rm m}^{-3}]\right)$")
    ax.set_ylabel(r"$\log_{10}\left(n_{Ne(+1)}\,[{\rm m}^{-3}]\right)$")
    ax.set_title(fr"$E_\parallel={ELECTRIC_FIELD}\,\rm V/m$")

    plt.tight_layout()
    plt.show()
