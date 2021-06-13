
import numpy as np
from scipy.optimize import least_squares

from . data import getADASIonizationData
from . helpers import averagedIonizationCrossSection


def fitKineticIonizationForSpecies(species, Z0, fittype):
    """
    Fits a kinetic ionization cross section for the given charge state Z0
    of the named ion species.

    :param species: Name of ion species.
    :param Z0:      Ion charge state.
    :param fittype: Method to use for fitting.
    """
    I_ADAS, Z, _, T, E_ion = getADASIonizationData(species)

    if Z0 >= Z:
        raise Exception("Invalid charge state specified: {}. Charge state must be strictly less than the atomic charge Z = {}.".format(Z0, Z))

    T_lower = 2
    T_upper = 100

    K1 = (T>T_lower).nonzero()[0][0]
    K2 = (T>T_upper).nonzero()[0]

    if K2.size == 0: K2 = T.size
    else: K2 = K2[0]

    I_A = I_ADAS[Z0,K1:K2,0]
    T_trunc = T[K1:K2]

    C_lo = 1e-4
    C_up = 200
    C_s  = 15
    betaStar_lo = 0
    betaStar_up = 1
    betaStar_s  = 0.5
    DI_lo = E_ion[0]
    DI_up = 100*DI_lo
    DI_s  = 10*DI_lo
    tolX  = 1e-6

    # Select fitting method
    outp = {'C1': None, 'C2': None, 'DI1': None, 'DI2': None, 'betaStar': None, 'beta2': None}
    if fittype == 'single':
        fo, goodness, output = _inner_fit(
            T_trunc.T, np.log(I_A), DI1=E_ion[Z0],
            fittype='single',
            lower=(C_lo, betaStar_lo), upper=(C_up, betaStar_up),
            guess=(C_s, betaStar_s), tolx=tolX
        )

        outp['C1'] = output[0]
        outp['DI1'] = E_ion[Z0]
        outp['betaStar'] = output[1]
    elif fittype == 'single_3p':
        fo, goodness, output = _inner_fit(
            T_trunc.T, np.log(I_A),
            fittype='single_3p',
            lower=(C_lo, betaStar_lo, E_ion[Z0]/3),
            upper=(C_up, betaStar_up, E_ion[Z0]*3),
            guess=(20, betaStar_s, E_ion[Z0]), tolx=tolX
        )

        outp['C1'] = output[0]
        outp['DI1'] = output[2]
        outp['betaStar'] = output[1]
    elif fittype == 'double':
        fo, goodness, output = _inner_fit(
            T_trunc.T, np.log(I_A), DI1=E_ion[Z0], beta2=0,
            fittype='double',
            lower=(C_lo, C_lo, betaStar_lo, DI_lo),
            upper=(C_up, C_up, betaStar_up, DI_up),
            guess=(C_s, C_s, betaStar_s, DI_s), tolx=tolX
        )
        outp['C1'] = output[0]
        outp['C2'] = output[1]
        outp['DI1'] = E_ion[Z0]
        outp['DI2'] = output[3]
        outp['betaStar'] = output[2]
        outp['beta2'] = 0
    else:
        raise Exception("Unrecognized 'fittype' specified: '{}'.".format(fittype))

    return fo, goodness, outp


def _inner_fit(x, y, lower, upper, guess, tolx, fittype='single', DI1=None, beta2=None):
    """
    Inner routine for fitting.
    """
    def f_single(x, T, y, DI1=None):
        C1 = x[0]
        betaStar = x[1]

        if DI1 is None:
            DI1 = x[2]

        return np.log(averagedIonizationCrossSection(T, C1, DI1, betaStar))-y

    def f_double(x, T, y, DI1=None, beta2=None):
        C1 = x[0]
        C2 = x[1]
        betaStar = x[2]
        DI2 = x[3]

        if DI1 is None:
            DI1 = x[4]
            beta2 = x[5]

        a1 = averagedIonizationCrossSection(T, C1, DI1, betaStar)
        a2 = averagedIonizationCrossSection(T, C2, DI2, beta2)
        return np.log(a1 + a2)-y

    if fittype in ['single', 'single_3p']:
        f = f_single
    elif fittype == 'double':
        f = f_double

    # Are DI1 and beta2 given explicitly, or are they free parameters?
    kwargs = {'T': x, 'y': y}
    if DI1 is not None: kwargs['DI1'] = DI1
    if beta2 is not None: kwargs['beta2'] = beta2

    result = least_squares(f, x0=guess, bounds=(lower, upper), ftol=0, xtol=tolx, kwargs=kwargs)

    return result, result.optimality, result.x


