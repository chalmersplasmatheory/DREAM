import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Ensure DREAM python package is importable when running from a repo checkout
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'py')))

from DREAM import DREAMSettings, runiface
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
import DREAM.Settings.Equations.BootstrapCurrent as BootstrapCurrent

# Timescales
tmax_TQ_max = 5e-3     # Maximum time for temporary thermal quench simulation (s)
tmax_CQ = 2e-1         # Total simulation time (s)
tmax = 50e-3

# Radial grid resolution
NR = 100

T0 = 12e3

tau_inv = 0

R0 = np.inf
a  = 1

def setStellaratorRadialGrid(ds, file, nr=NR, fb=1.1, load=False, visualize=False):
    """
    Set a stellarator magnetic field for the radial grid generator.

    Equilibrium provided by Stefan Buller, and can be found at: 
    https://sb0095.mycpanel.princeton.edu/QA/notes.html
    """
    ds.radialgrid.setType(RadialGrid.TYPE_STELLARATOR)

    if load: 
        ds.radialgrid.setStellarator(file, nr_equil=nr+1, ntheta_equil=65, nphi_equil=65, loadfilename=f'{file[:-3]}_settingsData')
    else: 
        ds.radialgrid.setStellarator(file, nr_equil=nr+1, ntheta_equil=65, nphi_equil=65, savefilename=f'{file[:-3]}_settingsData')

    ds.radialgrid.setNr(nr)
    ds.radialgrid.setMinorRadius(ds.radialgrid.num_stellarator.getMinorRadius())
    ds.radialgrid.setWallRadius(ds.radialgrid.num_stellarator.getMinorRadius()*fb)
    ds.radialgrid.setMajorRadius(ds.radialgrid.num_stellarator.getMajorRadius())

    global a, R0
    a = ds.radialgrid.num_stellarator.getMinorRadius()
    R0 = ds.radialgrid.num_stellarator.getMajorRadius()
    if visualize:
        ds.radialgrid.visualize(phi=np.linspace(0, 2 * np.pi / ds.radialgrid.nfp, 5)[:-1], nrho=5)


def getElectronTemperatureProfileEvolution(Tinit=T0, Tend=1, tau=0.5e-3):
    """  Returns the initial temperature profile. """
    r = np.linspace(0,1,101)
    T0_r = Tinit * (1 - 0.99 * r**2)
    t = np.linspace(0,tmax,1001)
    Te_mat = np.maximum(np.tile(T0_r - Tend * np.ones_like(T0_r), (1001, 1))*(np.exp(-t/tau).reshape((1001,1))) + Tend*np.ones_like(T0_r).reshape((1,101)), Tend)
    t = np.append(t,1)
    Te_mat = np.vstack((Te_mat, (Tend * np.ones_like(T0_r)).reshape((1,101))))
    return r*a, t, Te_mat
    

def getInitialElectronDensityProfile(n0=1e20, b=0.4):
    """  Returns the initial temperature profile. """
    r = np.linspace(0,1,101)
    c = (1 - 0.001**(1/b))
    profile = (1 - c * r**2)**b
    return r*a, n0*profile


def getInitialIonDensityProfile(ne0=1e20, b=0.4):
    """  Returns the initial temperature profile. """
    r, ne = getInitialElectronDensityProfile(ne0,b)

    return r*a, ne/2

def getInitialCurrentDensityProfile(ds, bootstrap=False, s=0, p=1, nr=NR):
    """ Returns the initial current density profile. """

    if bootstrap:
        r, tTe, Te = getElectronTemperatureProfileEvolution()
        ds.eqsys.T_cold.setType(Temperature.TYPE_PRESCRIBED)
        ds.eqsys.T_cold.setPrescribedData(Te[0,:], radius=r)

        ds.eqsys.E_field.setType(EField.TYPE_PRESCRIBED_OHMIC_CURRENT)
        ds.eqsys.j_ohm.setCurrentProfile(np.ones_like(r), radius=r, Ip0=10e6)
        ds.eqsys.j_bs.setMode(BootstrapCurrent.BOOTSTRAP_MODE_STELLARATOR)
        ds.eqsys.j_bs.setInitMode(BootstrapCurrent.BOOTSTRAP_INIT_MODE_TOTAL)

        ds.timestep.setTmax(1)
        ds.timestep.setNt(1)

        do_init = runiface(ds, quiet=True)
        j = do_init.eqsys.j_bs[0,:]*1.001
        r = do_init.grid.r[:]
        do_init.close()
    
    else:
        r = np.linspace(0,1,101)
        if s > 0:
            q = (1-s)/(1+s)
            j = np.flip(np.sin(np.pi*r**q)**p)
        else:
            q = (1+s)/(1-s)
            j = np.sin(np.pi*r**q)**p
        r *= a
    return r, j
