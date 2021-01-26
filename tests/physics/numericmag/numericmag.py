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


def constructMagneticField(Rp=2, Zp=0, a=0.5, nR=50, ntheta=51,
    rG=None, G=None, rpsi=None, psi=None,
    Delta=None, rDelta=None, kappa=None, rkappa=None,
    delta=None, rdelta=None):
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
    else: ipsi = UnivariateSpline(rpsi, psi, k=1, ext=3)

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

    return ipsi, theta, R, Z, Br, Bz, Bphi


