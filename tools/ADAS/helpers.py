# Implementation of various helper functions

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from scipy.integrate import quad
from scipy.special import kve
import sys


# Electron rest energy in eV
mc2 = 0.51099895000e6


def averagedIonizationCrossSection(T, C, DI_eV, betaStar):
    """
    """
    c = 299792458.0
    nT = T.size
    I_i = np.zeros(nT)

    def intg(p, temperature):
        pf = p**3 / np.sqrt(1+p**2)
        kic = kineticIonizationContribution(p, C, DI_eV, betaStar)
        fMe = maxwellJuttnerDistribution(p, 1, temperature)
        return pf*kic*fMe

    for k in range(nT):
        # SciPy's "quad()" seems to sometimes have problems with this integrand
        # over the infinite interval, so we split the integral into two parts:
        # one over an interval which should contain most, if not all, of the
        # interesting bits of the integrand, and one part covering the rest.
        pmax = np.sqrt((800*T[k]/mc2 + 1)**2 - 1)
        q = quad(lambda p : intg(p, T[k]), 0, pmax, epsabs=0)[0] + quad(lambda p : intg(p, T[k]), pmax, np.inf, epsabs=0)[0]
        I_i[k] = 4*pi*c*q

    return I_i


def maxwellJuttnerDistribution(p, n, T):
    """
    Evaluates a Maxwell-Jüttner distribution function at the
    given momentum, density and temperature, normalized such that
      
      integral( 4*pi*p^2 * fMe, 0, inf)
    
    yields the density ``n`` for all values of ``T``.
    """
    global mc2
    Theta = T / mc2

    tK2exp = 4*pi*Theta*kve(2, 1/Theta)

    gamma = np.sqrt(1+p**2)
    gMinus1 = p*p/(gamma+1)

    fMe = n/tK2exp * np.exp(-gMinus1/Theta)

    return fMe


def kineticIonizationContribution(p, C, DI_eV, betaStar):
    """
    Evaluates the  total electron impact ionization cross-section.

    :param p:        Incident electron momentum in units mc2.
    :param C:        Pre-factor (undetermined by the theory, or order ~1-10).
    :param DI_eV:    Ionization energy in eV.
    :param betaStar: Parameter which sets the near-threshold modification to
                     the cross-section.
    """
    global mc2
    a0  = 5.29e-11
    Ry  = 13.6 / mc2

    gamma = np.sqrt(1+p**2)
    Ek    = p*p/(gamma+1)
    DI    = DI_eV/mc2

    if DI <= 0:
        raise Exception('Invalid ionization energy provided (<=0)')
    
    U = Ek/DI

    if np.isscalar(U):  # handle vector momentum vector input
        if U > 1:
            I_nonRel = pi*a0**2*C*(Ry/DI)**2 * np.log(U)**(1+betaStar/U)/U
        else:
            I_nonRel = 0
    else:
        I_nonRel = pi*a0**2*C*(Ry/DI)**2 * np.log(U)**(1+betaStar/U)/U
        I_nonRel[np.where(U<=1)] = 0

    # v/c
    beta  = p/gamma
    # Fine structure constant
    alpha = 1/137
    # Expression appearing inside the log term of
    # the ultra-relativistic formula
    logArg    = p**2/(2*DI)
    if np.isscalar(U):
        if U > 1:
            I_rel = pi*a0**2*alpha**2 * C*(Ry/DI) * (np.log(logArg) - beta**2)
        else:
            I_rel = 0
    else:
        I_rel = pi*a0**2*alpha**2 * C*(Ry/DI) * (np.log(logArg) - beta**2)
        I_rel[np.where(U<=1)] = 0

    Ek_eV     = Ek*mc2
    S         = 1/(1+np.exp(1-Ek_eV*1e-5))
    I_kinetic = (1-S)*I_nonRel + S*I_rel

    return I_kinetic


if __name__ == '__main__':
    C = 3.024126380275995
    DI_eV = 13.598434490
    betaStar = 0.291350946928159
    T = 2.5113659405400823

    pmax = np.sqrt((800*T/mc2 + 1)**2 - 1)
    p = np.linspace(0, pmax, 1000)

    pf = p**3 / np.sqrt(1+p**2)
    kic, fMe = np.zeros(p.shape), np.zeros(p.shape)

    for i in range(p.size):
        kic[i] = kineticIonizationContribution(p[i], C, DI_eV, betaStar)
        fMe[i] = maxwellJuttnerDistribution(p[i], 1, T)

    print('II = [',end="")
    for i in range(p.size):
        print('{:.12e},'.format(pf[i]*kic[i]*fMe[i]), end="")

    print('];')
        
    plt.semilogy(p, pf*kic*fMe)
    plt.show()

