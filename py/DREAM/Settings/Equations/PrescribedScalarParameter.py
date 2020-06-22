# Helper class for prescribed parameters (i.e.
# radial profiles given at t=0)

import numpy as np

from . EquationException import EquationException


class PrescribedScalarParameter:
    

    def _setScalarData(self, data, times=0):
        """
        Set prescribed scalar data appropriately.
        """
        if np.isscalar(times):
            t = np.asarray([times])
        else: t = np.asarray(times)

        if np.isscalar(data):
            d = data*np.ones((t.size, ))
        else: d = np.asarray(data)

        return d, t


    def _verifySettingsPrescribedScalarData(self, name, data, times):
        """
        Verify the structure of the prescribed data.
        """
        if len(data.shape) != 1:
            raise EquationException("{}: Invalid number of dimensions in prescribed scalar data. Expected one dimension (times).".format(name))
        elif len(times.shape) != 1:
            raise EquationException("{}: Invalid number of dimensions in radial grid of prescribed scalar data. Expected one dimension.".format(name))
        elif data.shape[0] != times.size:
            raise EquationException("{}: Invalid size of prescribed data: {}. Expected {} elements."
                .format(name, data.shape[0], times.size))


