# Helper class for prescribed parameters (i.e.
# radial profiles given at t=0)

import numpy as np

from . EquationException import EquationException


class PrescribedInitialParameter:
    

    def _setInitialData(self, data, radius=0):
        """
        Set prescribed initial data appropriately.
        """
        if np.isscalar(radius):
            r = np.asarray([radius])
        else: r = np.asarray(radius)

        if np.isscalar(data):
            d = data*np.ones((r.size, ))
        else: d = np.asarray(data)

        return d, r


    def _verifySettingsPrescribedInitialData(self, name, data, radius):
        """
        Verify the structure of the prescribed data.
        """
        if len(data.shape) != 1:
            raise EquationException("{}: Invalid number of dimensions in prescribed initial data. Expected one dimension (radius).".format(name))
        elif len(radius.shape) != 1:
            raise EquationException("{}: Invalid number of dimensions in radial grid of prescribed initial data. Expected one dimension.".format(name))
        elif data.shape[0] != radius.size:
            raise EquationException("{}: Invalid size of prescribed data: {}. Expected {} elements."
                .format(name, data.shape[0], radius.size))


