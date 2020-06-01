# CODE CONDUCTIVITY
#
# Calculates the conductivity using the full linearized collision operator in
# DREAM and compares to the same calculation made with CODE.
#
############################################################################

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import scipy.constants
import sys

import DREAM
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMSettings import DREAMSettings

import DREAM.Settings.CollisionHandler as Collisions


def gensettings(T, Z=1, EEc=1e-2, n=5e19, pMax=5):
    """
    Generate appropriate DREAM settings.

    T:    Electron temperature.
    Z:    Effective charge of plasma.
    EEc:  Electric field (in units of critical electric field).
    n:    Electron density.
    pMax: Maximum momentum on computational grid.
    """
    c    = scipy.constants.c
    e    = scipy.constants.e
    eps0 = 8.85418782e-12
    me   = 9.10938e-31

    lnLambda = 14.9-0.5*np.log(n/1e20) + np.log(T/1e3)
    Ec = n*lnLambda*(e**3) / (4*np.pi*(eps0**2)*me*(c**2))

    ds = DREAMSettings()

    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

    ds.eqsys.E_field.setPrescribedData(EEc*Ec)
    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z)   # Imaginary ion with charge Z
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=n, rT0=0, T0=T)
    
    ds.hottailgrid.setNxi(10)
    ds.hottailgrid.setNp(500)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setNr(1)

    ds.timestep.setTmax(1e-4)
    ds.timestep.setNt(5)

    ds.save('input.h5')

    return ds


def loadCODE(filename):
    """
    Load conductivity values from CODE output.
    """
    T, Z, sigma = (None,)*3

    with h5py.File(filename, 'r') as f:
        T     = f['T'][:]
        Z     = f['Z'][:]
        sigma = f['sigma'][:]

    return T, Z, sigma


def runTZ(T, Z):
    """
    Run DREAM for the specified values of temperature and ion charge.
    """
    ds = gensettings(T=T, Z=Z)
    E  = ds.eqsys.E_field[0,0]
    print(E)

    do = DREAM.runiface(ds, 'output.h5')

    j = do.currentDensity(t=-1)
    sigma = j / E

    return sigma


def run():
    """
    Run the test.
    """
    workdir = pathlib.Path(__file__).parent.absolute()

    T, Z, CODEsigma = loadCODE('{}/CODE-conductivities.mat'.format(workdir))

    nZ, nT = T.shape
    sigma = np.zeros((nZ, nT))
    for i in range(0, nZ):
        for j in range(0, nT):
            sigma[i,j] = runTZ(T[i,j], Z[i,j])
    
    # Save
    with h5py.File('{}/DREAM-conductivities.h5'.format(workdir), 'w') as f:
        f['T'] = T
        f['Z'] = Z
        f['sigma'] = sigma

    # Compare conductivities
    legs = []
    legh = []
    hN = None
    for i in range(0, nT):
        clr = cmap(i/nT)

        h,  = plt.plot(Z[:,i], CODEsigma[:,i], color=clr, linewidth=2)
        hN, = plt.plot(Z[:,i], sigma[:,i], 'x', color=clr, markersize=10, markeredgewidth=3)

        legs.append(r'$T = {:.0f}\,\mathrm{{eV}}$'.format(T[0,i]))
        legh.append(h)

    legs.append('$\mathrm{DREAM}$')
    legh.append(hN)

    plt.xlabel(r'$Z$')
    plt.ylabel(r'$\sigma\ \mathrm{(S/m)}$')
    plt.legend(legh, legs)

    plt.show()



