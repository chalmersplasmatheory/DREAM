
import numpy as np

from . PlasmaParameters import getEc, getED, getTauEETh


def getConnorHastieRunawayRate(T, n, E, Zeff, with_corrections=True):
    """
    Calculates the runaway rate according to the formula given by Connor & Hastie.

    :param float T:               Electron temperature (eV).
    :param float n:               Electron density (m^-3).
    :param float E:               Electric field strength (V/m).
    :param float Zeff:            Plasma effective charge.
    :param bool with_corrections: If ``True``, includes the weak electric field limit corrections to the runaway rate.
    """
    if E == 0: return 0.0

    Ec = getEc(T, n)
    if Ec >= E: return 0.0

    ED = getED(T, n)
    tauEE = getTauEETh(T, n)

    EEc = E/Ec
    EED = E/ED

    # "Undetermined" factor (~1 is usually good according to simulations!)
    C = 1.0
    h, eta, lmbd = (1.0,)*3
    
    if with_corrections:
        h    = 1.0/(3*(EEc-1)) * (EEc + 2*(EEc-2)*np.sqrt(EEc/(EEc-1)) - (Zeff-7)/(Zeff+1))
        etaf = np.pi/2 - np.arcsin(1-2/EEc)
        eta  = EEc*EEc/(4*(EEc-1)) * etaf**2
        lmbd = 8*EEc**2 * (1-1/(2*EEc) - np.sqrt(1-1/EEc))

    alpha = -3/16*(1+Zeff)*h
    return C*n/tauEE * np.power(EED, alpha) * np.exp(-lmbd/(4*EED) - np.sqrt(eta*(1+Zeff)/EED))


