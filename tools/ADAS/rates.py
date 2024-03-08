# Load ADAS rate and create an spline object

import numpy as np
import scipy.interpolate
from .data import _loadADAS


def get_rate(species, rate_type, cache=True):
    """
    Load the ADAS rate for the specified species, of the given ADAS rate type.

    :param str species:  Name of ion species to load data for.
    :param str datatype: Name of datatype to load ('acd', 'scd', 'plt', 'prb', 'ccd')
    """
    req_shift = ['acd', 'ccd', 'prb']
    Z, n, T, x = _loadADAS(species=species, datatype=rate_type, cache=cache)

    return ADASRate(name=rate_type, x=x, n=n, T=T, Z=Z, species=species, shiftup=rate_type.lower() in req_shift)


class ADASRate:
    

    def __init__(self, name, x, n, T, Z, species, shiftup=False):
        """
        ADAS rate object.
        """
        self.Z = Z
        self.name = name
        self.species = species

        interp = []

        if shiftup:
            interp.append(None)

        for i in range(x.shape[0]):
            #interp.append(scipy.interpolate.interp2d(n, T, x[i,:], kind='linear', bounds_error=False))
            interp.append(scipy.interpolate.RectBivariateSpline(n, T, x[i,:].T, kx=3, ky=3))

        if len(interp) < Z+1:
            interp.append(None)

        self.interp = interp


    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)


    def eval(self, Z0, n, T):
        """
        Evaluate this rate for the given charge state, density and temperature.
        """
        if n <= 0:
            return 0

        r = self.interp[Z0]

        if r is None:
            raise ValueError(f"Cannot evaluate ADAS {self.name} rate for '{self.species}' with Z0 = {Z0}.")

        ln, lT = np.log10(n), np.log10(T)
        exp = r(ln, lT)

        if exp.size == 1:
            return 10**exp[0]
        else:
            return 10**exp


    def deriv_ne(self, Z0, n, T):
        """
        Evaluate the derivative of this rate w.r.t. the electron density.
        """
        return self._deriv(Z0=Z0, n=n, T=T, dn=1, dT=0)


    def deriv_Te(self, Z0, n, T):
        """
        Evaluate the derivative of this rate w.r.t. the electron temperature.
        """
        return self._deriv(Z0=Z0, n=n, T=T, dn=0, dT=1)


    def _deriv(self, Z0, n, T, dn=0, dT=0):
        """
        Evaluate a derivative of this rate.
        """
        if n <= 0:
            return 0

        if dn != 0 and dT != 0:
            raise Exception("Cannot evaluate mixed derivatives of rates.")
        elif dn not in [0,1] and dT not in [0,1]:
            raise Exception("Cannot evaluate higher than first-order derivatives of rates.")

        r = self.interp[Z0]

        if r is None:
            raise ValueError(f"Cannot evaluate ADAS {self.name} rate for '{self.species}' with Z0 = {Z0}.")

        ln, lT = np.log10(n), np.log10(T)
        dI = r.partial_derivative(dn, dT)
        exp = r(ln, lT)
        expd = dI(ln, lT)

        if dn == 1:
            iv = 1/n
        else:
            iv = 1/T

        if exp.size == 1:
            return 10**exp[0] * expd * iv
        else:
            return 10**exp * expd * iv

            
