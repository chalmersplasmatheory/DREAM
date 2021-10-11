#!/usr/bin/env python3
#
# TOKAMAK:   ASDEX Upgrade
################################################################

import numpy as np
import scipy.constants

try:
    import DREAM.Settings.RadialGrid as RadialGrid
except:
    import sys
    sys.path.append('/run/media/mathias/BlackTools/Fusion/software/DREAM/py')
    sys.path.append('/home/ola/svn/runaway/ola/DREAM/py')


    import DREAM.Settings.RadialGrid as RadialGrid

import visualizeTokamak


# Plasma parameters
a  = 0.5        # minor radius (m)
b  = 0.55       # wall radius (m)
B0 = 2.5        # on-axis magnetic field strength (T)
R0 = 1.65       # major radius (m)
Ip = 800e3      # target plasma current (A)
j0 = 1

ne0 = 2.6e19    # central electron density (m^-3)
Te0 = 5.8e3     # central electron temperature (eV)


def setMagneticField(ds, nr=40, visualize=False, rGridNonuniformity=1, toroidal=True):
    """
    Set the radial grid of the given DREAMSettings object to correspond to a
    typical magnetic field in ASDEX Upgrade.
    """
    global R0, a, b, B0, Ip

    # Input radial grid resolution
    NR = 101

    # Radial grid for analytical magnetic field
    r = np.linspace(0, a, NR)

    #if toroidal:
    # Elongation profile
    kappa = np.linspace(1.0, 1.15, NR)

    # Poloidal flux radial grid
    mu0 = scipy.constants.mu_0
    psi_p = -mu0 * Ip * (1-(r/a)**2) * a

    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
    ds.radialgrid.setWallRadius(b)

    if toroidal:
        ds.radialgrid.setMajorRadius(R0)
    else:
        ds.radialgrid.setMajorRadius(np.inf)
        psi_p = np.zeros(r.shape)

    q = rGridNonuniformity
    r_f = np.linspace(0, a**q,nr+1)**(1/q)
    r_f = r_f * a/r_f[-1] # correct for roundoff
    ds.radialgrid.setCustomGridPoints(r_f)

    ds.radialgrid.setShaping(psi=psi_p, rpsi=r, GOverR0=B0, kappa=kappa, rkappa=r)

    if visualize:
        ds.radialgrid.visualize(ntheta=200)
    """
    else:   # Cylindrical
        ds.radialgrid.setType(RadialGrid.TYPE_CYLINDRICAL)
        ds.radialgrid.setB0(B0)
        ds.radialgrid.setMinorRadius(a)
        ds.radialgrid.setWallRadius(b)
        ds.radialgrid.setNr(nr)
    """


def getInitialTemperature(r=None, nr=100):
    """
    Returns the initial temperature profile.

    :param r:      Minor radius (m). May be either scalar or numpy.ndarray.
    :param int nr: Radial grid resolution to use if ``r = None``.
    """
    global a, Te0

    if r is None:
        r = np.linspace(0, a, nr)

    #return r, Te0 * (1 - 0.99*(r/a)**2)
    c1 = 0.4653
    return r, Te0 * np.exp(-((r/a)/c1)**2)


def getInitialDensity(r=None, nr=100):
    """
    Returns the initial density profile.

    :param r:      Minor radius (m). May be either scalar or numpy.ndarray.
    :param int nr: Radial grid resolution to use if ``r = None``.
    """
    global a, ne0

    if r is None:
        r = np.linspace(0, a, nr)

    #return r, ne0 * (1 - 0.99*(r/a)**2)
    return r, ne0 * np.polyval([-0.1542, -0.01115, 1.], r/a)


def getCurrentDensity(r=None, nr=100):
    """
    Returns the initial current density profile.

    :param r:      Minor radius (m). May be either scalar or numpy.ndarray.
    :param int nr: Radial grid resolution to use if ``r = None``.
    """
    global a

    if r is None:
        r = np.linspace(0, a, nr)

    j0 = 1

    #return r, j0 * (1 - (1-0.001**(1/0.41))*(r/a)**2)**0.41
    return r, j0 * (1 - (r/a)**4)**1.5


if __name__ == '__main__':
    visualizeTokamak.visualize(__import__(__name__), nticks=[0, 1, 2], Tticks=[0, 2, 4, 6])  #sys.modules[__name__])


