# Helper class for prescribed parameters

import numpy as np

from . EquationException import EquationException


class PrescribedParameter:
    

    def _setPrescribedData(self, data, radius=0, times=0):
        """
        Set prescribed data appropriately.
        """
        if np.isscalar(radius):
            r = np.asarray([radius])
        else: r = np.asarray(radius)

        if np.isscalar(times):
            t = np.asarray([times])
        else: t = np.asarray(times)

        if np.isscalar(data):
            d = data*np.ones((t.size, r.size))
        else: d = np.asarray(data)

        return d, r, t


    def _verifySettingsPrescribedData(self, name, data, radius, times):
        """
        Verify the structure of the prescribed data.
        """
        if len(data.shape) != 2:
            raise EquationException("{}: Invalid number of dimensions in prescribed data. Expected 2 dimensions (time x radius).".format(name))
        elif len(times.shape) != 1:
            raise EquationException("{}: Invalid number of dimensions in time grid of prescribed data. Expected one dimension.".format(name))
        elif len(radius.shape) != 1:
            raise EquationException("{}: Invalid number of dimensions in radial grid of prescribed data. Expected one dimension.".format(name))
        elif data.shape[0] != times.size or data.shape[1] != radius.size:
            raise EquationException("{}: Invalid dimensions of prescribed data: {}x{}. Expected {}x{} (time x radius)."
                .format(name, data.shape[0], data.shape[1], times.size, radius.size))


