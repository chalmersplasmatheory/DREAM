#!/usr/bin/env python3
#
# This example generates a tokamak magnetic field using the analytica toroidal
# model, implemented separately in DREAM, and stores it as a numeric magnetic
# field. Two simulations are then run, utilizing both the analytic and
# corresponding numeric magnetic fields, and the results are compared.
#

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import scipy.constants
from scipy.interpolate import UnivariateSpline
import sys

import dreamtests

import DREAM
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMSettings import DREAMSettings
import DREAM.GeriMap as GeriMap

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.RunawayElectrons as Runaways


Rp = 2.0        # Major radius
a  = 0.55       # Minor radius
B0 = 6.2


def getShapeProfiles(nr=20):
    """
    Returns the shape profiles for magnetic field to use.
    """
    global a, B0, Rp

    r = np.linspace(0, a, nr)

    rG, G         = r, B0 * np.ones(r.shape)
    rDelta, Delta = r, np.linspace(0, 0.05*a, nr)
    rkappa, kappa = r, np.linspace(1, 1.4, nr)
    rdelta, delta = r, np.linspace(0, 0.05, nr)

    IpRef = 7e6
    mu0   = scipy.constants.mu_0
    rpsi  = r
    psi   = -mu0 * IpRef * (1-(rpsi/a)**2) * a

    return rG, G, rpsi, psi, rDelta, Delta, rkappa, kappa, rDelta, Delta


def constructMagneticField(Rp=2, Zp=0, a=0.5, nR=50, ntheta=51,
    rG=None, G=None, rpsi=None, psi=None,
    Delta=None, rDelta=None, kappa=None, rkappa=None,
    delta=None, rdelta=None, retdict=False):
    """
    Construct the numeric magnetic field.
    """
    r = np.linspace(0, a, nR+1)[1:]
    theta = np.linspace(0, 2*np.pi, ntheta+1)[:-1]

    mgR, mgT = np.meshgrid(r, theta)

    iG, ipsi, iDelta, ikappa, idelta = (None,)*5

    if G is None: raise Exception('The toroidal magnetic field function must be specified.')
    else: iG = UnivariateSpline(rG, G, k=1, ext=3)

    if psi is None: raise Exception('The poloidal flux must be specified.')
    else: ipsi = UnivariateSpline(rpsi, psi/(2*np.pi), k=1, ext=3)

    if Delta is None: iDelta = UnivariateSpline([0, 1], [0, 0], k=1, ext=3)
    else: iDelta = UnivariateSpline(rDelta, Delta, k=1, ext=3)

    if kappa is None: ikappa = interp1d([0], [1], k=1, ext=3)
    else: ikappa = UnivariateSpline(rkappa, kappa, k=1, ext=3)

    if delta is None: idelta = interp1d([0], [0], k=1, ext=3)
    else: idelta = UnivariateSpline(rdelta, delta, k=1, ext=3)

    R = iDelta(r) + r*np.cos(theta + idelta(r)*np.sin(theta))
    Z = r*ikappa(r)*np.sin(theta)

    gradPsi = ipsi.derivative()

    dRdr, dZdr = np.zeros(R.shape), np.zeros(Z.shape)
    dRdr[:,:-1] = np.diff(R) / np.diff(r)
    dZdr[:,:-1] = np.diff(Z) / np.diff(r)

    dRdr[:,-1] = dRdr[:,-2]
    dZdr[:,-1] = dZdr[:,-2]

    Bphi = iG(r)
    Br   = -gradPsi(r) / R * dZdr / np.sqrt(dRdr**2 + dZdr**2)
    Bz   =  gradPsi(r) / R * dRdr / np.sqrt(dRdr**2 + dZdr**2)

    if retdict:
        return {'Rp': Rp, 'Zp': Zp, 'psi': ipsi(r), 'theta': theta, 'R': R, 'Z': Z, 'Br': Br, 'Bz': Bz, 'Bphi': Bphi}
    else:
        return Rp, Zp, ipsi(r), theta, R, Z, Br, Bz, Bphi


def generateSettings(analyticB=False):
    """
    Generates a DREAMSettings object for the test

    :param bool analyticB: If ``True``, uses an analytic magnetic field. Otherwise a numeric magnetic field is used.
    """
    T = 3e3     # eV
    E = 2       # V/m
    n = 5e19    # m^-3
    yMax = 20   # thermal momentum

    betaTh = DREAM.Formulas.getNormalizedThermalSpeed(T)
    pMax = yMax * betaTh
    Ec = DREAM.Formulas.getEc(T, n)

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_THERMAL

    ds.eqsys.E_field.setPrescribedData(E)
    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   # Imaginary ion with charge Z
    ds.eqsys.n_cold.setPrescribedData(n)
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=n, rT0=0, T0=T)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_QUICK)
    
    ds.hottailgrid.setNxi(20)
    ds.hottailgrid.setNp(100)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    # Get magnetic field shaping parameters
    rG, G, rpsi, psi, rDelta, Delta, rkappa, kappa, rDelta, Delta = getShapeProfiles()
    if analyticB:
        ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
        ds.radialgrid.setShaping(psi=psi, rpsi=rpsi, G=G, rG=rG, kappa=kappa, rkappa, delta=delta, rdelta=rdelta, Delta=Delta, rDelta)
    else:
        raise Exception("Support for numeric magnetic fields not yet implemented.")

    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a)
    ds.radialgrid.setMajorRadius(Rp)
    ds.radialgrid.setNr(10)

    tMax0 = pMax*Ec / E
    ds.timestep.setTmax(.9*tMax0)
    ds.timestep.setNt(nTimeSteps)

    ds.other.include('fluid/runawayRate', 'fluid/gammaDreicer')

    return ds


def run(args):
    """
    Run the test.
    """
    success = True

    return success

