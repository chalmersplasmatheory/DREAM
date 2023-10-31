#!/usr/bin/env python3
# Definitions of ITER parameters


import numpy as np

try:
    import DREAM.Settings.RadialGrid as RadialGrid
except:
    import sys
    sys.path.append('/run/media/mathias/BlackTools/Fusion/software/DREAM/py')

    import DREAM.Settings.RadialGrid as RadialGrid

import visualizeTokamak


# ITER parameters
a  = 1.79    # minor radius (m)
b  = 2.00    # minor radius of tokamak wall (m)
B0 = 5.3     # toroidal magnetic field on-axis (T)
R0 = 6.2     # major radius (m)
Ip = 15e6    # Target plasma current (A)
j0 = 1.69e6  # current density (A/m^2) (yielding a plasma current Ip=15 MA in a circular plasma)

# Simulation parameters
ne0 = 1e20  # electron density (m^-3)
t0  = 1e-3  # xponential decay time (s)


def setMagneticField(ds, nr=40, visualize=False):
    """
    Set an ITER-like magnetic field for the AnalyticB radial grid generator.
    """
    global R0, a, b, B0

    # Radial grid resolution
    NR = 101

    # Radial grid for analytical magnetic field
    r = np.linspace(0, a, NR)
    
    # Elongation profile (fit)
    kappa = np.polyval([0.13408883,-0.35950981,0.37803907,1.17567818], r)

    # Triangularity profile
    delta = 0.35*(r/a)

    # Poloidal flux radial grid
    psi_r = [0.,0.15729669,0.21956937,0.26821122,0.30894407,0.34456487,0.37669775,0.40612651,0.43344149,0.45903068,0.48318347,0.50611202,0.52798629,0.54893893,0.56907693,0.58848811,0.60724724,0.62541499,0.64304163,0.66016932,0.67683145,0.69305626,0.70886366,0.72427377,0.73930629,0.7539885,0.76834598,0.78240391,0.7961856,0.80971195,0.82300201,0.83607276,0.84893872,0.86161257,0.87410597,0.88642996,0.89859466,0.91060962,0.92248405,0.93422631,0.94584417,0.95734403,0.96873219,0.98001427,0.99119554,1.00228094,1.01327516,1.02418297,1.03500895,1.04575739,1.05643194,1.06703588,1.07757233,1.08804434,1.09845488,1.10880682,1.11910278,1.129345,1.13953555,1.14967642,1.15976976,1.16981777,1.17982258,1.18978613,1.19971007,1.20959591,1.21944516,1.22925925,1.23903956,1.24878743,1.25850418,1.26819105,1.27784926,1.28747997,1.29708427,1.30666306,1.31621725,1.32574772,1.33525552,1.34474181,1.3542077,1.36365421,1.37308209,1.38249198,1.39188451,1.40126024,1.41061957,1.41996291,1.42929063,1.43860331,1.44790164,1.45718634,1.46645808,1.47571723,1.48496403,1.49419869,1.50342147,1.51263288,1.52183344,1.53102371,1.54020416,1.54937504,1.55853657,1.56768898,1.57683252,1.58596753,1.59509438,1.60421339,1.61332461,1.62242774,1.63152244,1.64060842,1.64968555,1.65875389,1.6678135,1.67686424,1.6859037,1.69492821,1.7039341,1.71291765,1.72187424,1.73079888,1.73968684,1.74853379,1.75733735,1.76609582,1.77480767,1.78347367,1.79100379]
    # Poloidal flux
    psi_p = [0.,0.02707751,0.05415503,0.08123254,0.10831006,0.13538757,0.16246508,0.1895426,0.21662011,0.24369763,0.27077514,0.29785265,0.32493017,0.35200768,0.3790852,0.40616271,0.43324022,0.46031774,0.48739525,0.51447276,0.54155028,0.56862779,0.59570531,0.62278282,0.64986033,0.67693785,0.70401536,0.73109288,0.75817039,0.7852479,0.81232542,0.83940293,0.86648045,0.89355796,0.92063547,0.94771299,0.9747905,1.00186802,1.02894553,1.05602304,1.08310056,1.11017807,1.13725559,1.1643331,1.19141061,1.21848813,1.24556564,1.27264315,1.29972067,1.32679818,1.3538757,1.38095321,1.40803072,1.43510824,1.46218575,1.48926327,1.51634078,1.54341829,1.57049581,1.59757332,1.62465084,1.65172835,1.67880586,1.70588338,1.73296089,1.76003841,1.78711592,1.81419343,1.84127095,1.86834846,1.89542598,1.92250349,1.949581,1.97665852,2.00373603,2.03081354,2.05789106,2.08496857,2.11204609,2.1391236,2.16620111,2.19327863,2.22035614,2.24743366,2.27451117,2.30158868,2.3286662,2.35574371,2.38282123,2.40989874,2.43697625,2.46405377,2.49113128,2.5182088,2.54528631,2.57236382,2.59944134,2.62651885,2.65359637,2.68067388,2.70775139,2.73482891,2.76190642,2.78898393,2.81606145,2.84313896,2.87021648,2.89729399,2.9243715,2.95144902,2.97852653,3.00560405,3.03268156,3.05975907,3.08683659,3.1139141,3.14099162,3.16806913,3.19514664,3.22222416,3.24930167,3.27637919,3.3034567,3.33053421,3.35761173,3.38468924,3.41176676,3.43884427,3.46592178]

    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
    ds.radialgrid.setMajorRadius(R0)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(b)
    ds.radialgrid.setNr(nr)

    ds.radialgrid.setShaping(psi=psi_p, rpsi=psi_r, G=B0, kappa=kappa, rkappa=r, delta=delta, rdelta=r)
    
    if visualize:
        ds.radialgrid.visualize(ntheta=200)


def getInitialTemperature(r=None, nr=100):
    """
    Returns the initial temperature profile.

    :param r:      Minor radius (m). May be either scalar or numpy.ndarray.
    :param int nr: Radial grid resolution to use if ``r = None``.
    """
    global a

    if r is None:
        r = np.linspace(0, a, nr)

    # Parabolic temperature profile
    #return r, 20e3 * (1 - (r/a))**2
    return r, 20e3 * (1 - 0.99*(r/a)**2)


def getFinalTemperature(r=None, nr=100):
    """
    Returns the final temperature profile.

    :param r:      Minor radius (m). May be either scalar or numpy.ndarray.
    :param int nr: Radial grid resolution to use if ``r = None``.
    """
    global a

    if r is None:
        r = np.linspace(0, a, nr)
    elif type(r) == float or type(r) == int:
        r = np.array([float(r)])

    # Flat temperature profile
    return r, 50*np.ones(r.shape)


def getTemperatureEvolution(tau0, nt=100):
    """
    Returns the spatiotemporal temperature profile

    :param float tau0: Time decay parameter.
    """
    r, T0 = getInitialTemperature()
    _, Tf = getFinalTemperature(r)

    # Let the temperature drop until the central temperature is T~100 eV
    #maxt = -tau0 * np.log((100 - Tf[0]) / (T0[0] - Tf[0]))
    maxt = 6e-3

    # Generate time grid
    t = np.linspace(0, maxt, nt).reshape((nt,1))

    # Reshape
    nr = r.size
    T0 = T0.reshape((1,nr))
    Tf = Tf.reshape((1,nr))

    # Calculate temperature evolution
    T = Tf + (T0 - Tf) * np.exp(-t/tau0)
    t = t.reshape((nt,))

    return t, r, T


def getInitialDensity(r=None, nr=100):
    """
    Returns the initial electron density profile.
    """
    global a, ne0

    if r is None:
        r = np.linspace(0, a, nr)

    return r, ne0 * np.ones(r.shape)


def getCurrentDensity(r=None, nr=100):
    """
    Returns the initial current density profile.

    :param r:      Minor radius (m). May be either scalar or numpy.ndarray.
    :param int nr: Radial grid resolution to use if ``r = None``.
    """
    global a, j0

    if r is None:
        r = np.linspace(0, a, nr)

    #return r, j0 * (1 - (r/a)**0.41)
    return r, j0 * (1 - (1-0.001**(1/0.41))*(r/a)**2)**0.41
    

if __name__ == '__main__':
    visualizeTokamak.visualize(__import__(__name__))  #sys.modules[__name__])


