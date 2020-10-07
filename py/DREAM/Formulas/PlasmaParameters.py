
import numpy as np
import scipy.constants


def getCoulombLogarithm(T, n):
    """
    Calculates the Coulomb logarithm according to the formula given in
    Wesson's book "Tokamaks".

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    return 14.9 - 0.5*np.log(n / 1e20) + np.log(T / 1e3)


def getEc(T, n):
    """
    Calculates the Connor-Hastie critical electric field, below which no
    runaway electrons can be generated.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    logLambda = getCoulombLogarithm(T, n)

    c = scipy.constants.c
    e = scipy.constants.e
    me = scipy.constants.m_e
    eps0 = scipy.constants.epsilon_0

    return (n*logLambda*e**3) / (4*np.pi*eps0**2 * me * c**2)


def getConnorHastieCriticalField(T, n): return getEc(T, n)


def getED(T, n):
    """
    Calculates the Dreicer electric field at the given plasma temperature and
    density, giving the electric field at which all electrons are accelerated
    to the runaway region.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    c = scipy.constants.c
    e = scipy.constants.e
    me = scipy.constants.m_e

    Ec = getEc(T, n)

    return Ec * me * c**2 / (e*T)


def getDreicerElectricField(T, n): return getED(T, N)


def getTauEETh(T, n):
    """
    Calculates the thermal electron-electron collision frequency.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    mc2  = scipy.constants.physical_constants['electron mass energy equivalent in MeV'][0] * 1e6
    vth2 = 2*T/mc2

    return getTauEERel(T, n) * vth2*np.sqrt(vth2)


def getThermalElectronCollisionFrequency(T, n): return getTauEETh(T, n)


def getTauEERel(T, n):
    """
    Calculates the relativistic electron-electron collision frequency.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    c = scipy.constants.c
    e = scipy.constants.e
    me = scipy.constants.m_e
    r0 = scipy.constants.physical_constants['classical electron radius'][0]

    C = 4*np.pi * r0**2 * c
    logLambda = getCoulombLogarithm(T, n)

    return 1/(logLambda * n * C)


def getRelativisticElectronCollisionFrequency(T, n): return getTauEERel(T, n)


