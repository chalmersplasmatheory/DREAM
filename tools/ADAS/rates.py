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

    return ADASRate(x=x, n=n, T=T, Z=Z, species=species, shiftup=rate_type.lower() in req_shift)


class ADASRate:
    

    def __init__(self, x, n, T, Z, species, shiftup=False):
        """
        ADAS rate object.
        """
        self.Z = Z
        self.species = species

        interp = []

        if shiftup:
            interp.append(None)

        for i in range(x.shape[0]):
            interp.append(scipy.interpolate.interp2d(n, T, x[i,:], kind='linear', bounds_error=False))

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

            
