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
    

