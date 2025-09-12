import numpy as np
import sys

try:
    import DREAM
except ModuleNotFoundError:
    sys.path.append('~/Skrivbord/DREAM/py')
    import DREAM

import DREAM.Settings.RadialGrid as RadialGrid
from DREAM import DREAMOutput

# Need major radius, minor radius, magnetic energy, magnetic field
# Temperature, electron density
# What is the scenario they are most worried about
#   - Magnet quench?
#   - Radiative collapse due to impurities
#   - MHD instabillities? Probably not, should be optimized away, but maybe close to MHD limit?
# Have you done anything on this so far
#

# ITER parameters
R0 = 5.5             # major radius (m)
a  = 0.5             # minor radius (m)
b  = 0.6             # minor radius of tokamak wall (m)
B0 = 3               # toroidal magnetic field on-axis (T)
Ip = 10e6            # Target plasma current (A)

# Timescales
tmax_TQ_max = 2e-3     # Maximum time for temporary thermal quench simulation (s)
tmax_CQ = 20e-3         # Total simulation time (s)

# Simulation parameters
ne0 = 1e20  # electron density (m^-3)
T_initial = 5e3
f_T0_o = 0.99

# Radial grid resolution
NR = 101

def setRadialGrid(ds, nr=NR):
    """
    Set an ITER-like magnetic field for the radial grid generator.

    Source: István Pusztai.
    """
    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)

    ds.radialgrid.setMajorRadius(R0)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(b)
    ds.radialgrid.setNr(nr)
    #ds.radialgrid.setB0(B0)
    r = np.linspace(0, a, nr)
    ds.radialgrid.setShaping(psi=r, rpsi=r**2, GOverR0=B0, kappa=1, Delta=0, delta=0)
    #ds.radialgrid.visualize(ntheta=200, show=True)
    '''
    # Radial grid for analytical magnetic field
    r = np.linspace(0, a, nr)

    def f(a0, a1, a2, a3, a4):
        return np.polyval(np.flip([a0, a1, a2, a3, a4]), r)
    
    # Shaping parameters
    ds.radialgrid.setShapeParameter('kappa', r=r, data=f(1.5, 0, 0, 0, .02))
    ds.radialgrid.setShapeParameter('delta', r=r, data=f(0, 0, .035, 0, .017))
    ds.radialgrid.setShapeParameter('Delta', r=r, data=f(0, 0, -.00658678, -.00634124, 0))

    # Poloidal magnetic flux
    ds.radialgrid.setShapeParameter('psi_p0', r=r, data=f(0, 0, .794102, -0.117139, 0))

    # Toroidal magnetic field function
    ds.radialgrid.setShapeParameter('GOverR0', r=r, data=f(B0, 0, -.310741, .107719, 0))
    #'''

def getInitialTemperatureProfile(nr=NR):
    """  Returns the initial temperature profile. """
    r = np.linspace(0, a, nr)
    T = T_initial * (1 - f_T0_o * (r/a)**2)
    return r, T

def getInitialDensityProfile(nr=NR):
    """  Returns the initial temperature profile. """
    r = np.linspace(0, a, nr)
    n = ne0*(1 - f_T0_o * (r/a)**2)#np.ones(nr)
    return r, n

def getInitialCurrentDensityProfile(nr=NR):
    """ Returns the initial current density profile. """
    r = np.linspace(0, a, nr)
    j2 = 0.41
    j1 = (1 - 0.001**(1./j2))
    j = (1 - j1*(r/a)**2)**j2
    do = DREAMOutput('output_BS.h5')
    return do.grid.r[:], do.eqsys.j_bs[0,:] + j*1e-6