import numpy as np
import scipy.constants 
import scipy.special

def getRelativisticMaxwellJuttnerDistribution(p,T,n):
    """
    Calculates the relativistic Maxwell-Jüttner distribution.

    :param float p: Momentum in units of m_e*c
    :param float T: Plasma temperature (eV)
    :param float n: Plasma density (m^-3)
    """
    mc2inEV = scipy.constants.m_e*scipy.constants.c*scipy.constants.c / scipy.constants.elementary_charge
    Theta = T / mc2inEV
    K2scaled = scipy.special.kve(2,1/Theta) # scaled bessel function 
    tK2exp = 4*np.pi*Theta*K2scaled # normalization factor
    p2 = p**2
    g = np.sqrt(1+p2)
    gMinus1 = p2/(g+1) 
    expTerm = np.exp(-gMinus1/Theta)

    return n * expTerm / tK2exp


def getAvalancheDistribution(p, xi, E, Z, nre=1, logLambda=15):
    """
    Evaluates the analytical avalanche distribution function according to
    equation (2.17) of [Embreus et al, J. Plasma Phys. 84 (2018)].

    :param p:         Momentum grid on which to evaluate the distribution.
    :param xi:        Pitch grid on which to evaluate the distribution.
    :param E:         Electric field strength (normalized to the Connor-Hastie field, Ec).
    :param Z:         Plasma total charge (= 1/ne_tot * sum_i ni * Zi^2)
    :param nre:       Runaway electron density.
    :param logLambda: Coulomb logarithm.
    """
    if p.ndim == 1:
        P, XI = np.meshgrid(p, xi)
    else:
        P, XI = p, xi

    c = scipy.constants.c
    m_e = scipy.constants.m_e

    g = np.sqrt(1+P**2)
    A = (E+1) / (Z+1) * g
    cZ = np.sqrt(5+Z)
    g0 = cZ*logLambda

    pf = m_e*c * nre * A / (2*np.pi*m_e*c*g0*P**2) / (1-np.exp(-2*A))
    f = pf * np.exp(-g/g0 - A*(1-XI))

    return f


