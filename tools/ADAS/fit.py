
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
import warnings

from . data import getIonizationData
from . helpers import averagedIonizationCrossSection


# Disable useless warning from 'scipy.optimize.least_squares()'
warnings.filterwarnings("ignore", message="Setting `ftol` below")
warnings.filterwarnings("ignore", message="divide by zero encountered in log")


def fitKineticIonizationForSpecies(species, Z0, fittype, T_lower=2, T_upper=100):
    """
    Fits a kinetic ionization cross section for the given charge state Z0
    of the named ion species.

    :param species: Name of ion species.
    :param Z0:      Ion charge state.
    :param fittype: Method to use for fitting.
    """
    I_ADAS, Z, _, T, E_ion = getIonizationData(species)

    if Z0 >= Z:
        raise Exception("Invalid charge state specified: {}. Charge state must be strictly less than the atomic charge Z = {}.".format(Z0, Z))

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


def evaluateAveragedCrossSection(T, C1, DI1, betaStar, method='single', C2=None, DI2=None, beta2=None):
    """
    Evaluate the averaged ionization cross section using the specified
    evaluation method and provided fit parameters.
    """
    if method in ['single', 'single_3p']:
        return averagedIonizationCrossSection(T, C1, DI1, betaStar)
    elif method == 'double':
        a1 = averagedIonizationCrossSection(T, C1, DI1, betaStar)
        a2 = averagedIonizationCrossSection(T, C2, DI2, beta2)
        return a1 + a2
    else:
        raise Exception("Unrecognized fitting method '{}'.".format(method))


def _inner_fit(x, y, lower, upper, guess, tolx, fittype='single', DI1=None, beta2=None):
    """
    Inner routine for fitting.
    """
    def f_single(x, T, y, DI1=None):
        C1 = x[0]
        betaStar = x[1]

        if DI1 is None:
            DI1 = x[2]

        v = np.log(evaluateAveragedCrossSection(T=T, C1=C1, DI1=DI1, betaStar=betaStar))-y

        # Remove points with ICS identically zero
        v[np.where(np.isinf(v))] = 0

        return v

    def f_double(x, T, y, DI1=None, beta2=None):
        C1 = x[0]
        C2 = x[1]
        betaStar = x[2]
        DI2 = x[3]

        if DI1 is None:
            DI1 = x[4]
            beta2 = x[5]

        #a1 = averagedIonizationCrossSection(T, C1, DI1, betaStar)
        #a2 = averagedIonizationCrossSection(T, C2, DI2, beta2)
        v = np.log(evaluateAveragedCrossSection(T=T, C1=C1, C2=C2, DI1=DI1, DI2=DI2, betaStar=betaStar, beta2=beta2, method='double'))-y

        # Remove points with ICS identically zero
        v[np.where(np.isinf(v))] = 0

        return v

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


