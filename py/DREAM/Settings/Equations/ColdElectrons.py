
import matplotlib.pyplot as plt
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . PrescribedScalarParameter import PrescribedScalarParameter
from . UnknownQuantity import UnknownQuantity


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2


class ColdElectrons(UnknownQuantity,PrescribedParameter, PrescribedScalarParameter):
    
    def __init__(self, settings, ttype=TYPE_SELFCONSISTENT, density=None, radius=None, times=None):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.setType(ttype=ttype)

        self.density = None
        self.radius  = None
        self.times   = None

        if (ttype == TYPE_PRESCRIBED) and (density is not None) and (radius is not None) and (times is not None):
            self.setPrescribedData(density=density, radius=radius, times=times)
        

    ###################
    # SETTERS
    ###################
    def setPrescribedData(self, density, radius=0, times=0):
        """
        Prescribe the cold electron density in time and space.

        :param density: Cold electron density (2D array (nt, nr) or scalar (=> constant and uniform in time and space))
        :param radius:  Radial grid on which the cold electron density is defined.
        :param times:   Time grid on which the cold electron density is defined.
        """
        _data, _rad, _tim = self._setPrescribedData(density, radius, times)
        self.density = _data
        self.radius = _rad
        self.times  = _tim

        self.setType(TYPE_PRESCRIBED)

        self.verifySettingsPrescribedData()


    def setType(self, ttype):
        """
        Sets the type of equation to use for evolving the cold electron density.

        :param int ttype: Flag indicating how to evolve the cold electron density.
        """
        if ttype == TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == TYPE_SELFCONSISTENT:
            self.type = ttype
        else:
            raise EquationException("n_cold: Unrecognized cold electron density type: {}".format(self.type))


    def fromdict(self, data):
        """
        Set all options from a dictionary.

        :param dict data: List of settings to load.
        """
        self.type = data['type']

        if self.type == TYPE_PRESCRIBED:
            self.density = data['data']['x']
            self.radius  = data['data']['r']
            self.times   = data['data']['t']
        elif self.type == TYPE_SELFCONSISTENT:
            pass
        else:
            raise EquationException("n_cold: Unrecognized cold electron density type: {}".format(self.type))

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }

        if self.type == TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.density,
                'r': self.radius,
                't': self.times
            }
        elif self.type == TYPE_SELFCONSISTENT:
            pass
        else:
            raise EquationException("n_cold: Unrecognized cold electron density type: {}".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.type == TYPE_PRESCRIBED:
            if type(self.density) != np.ndarray:
                raise EquationException("n_cold: Density prescribed, but no density data provided.")
            elif type(self.times) != np.ndarray:
                raise EquationException("n_cold: Density prescribed, but no time data provided, or provided in an invalid format.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("n_cold: Density prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedData()
        elif self.type == TYPE_SELFCONSISTENT:
            # Nothing todo
            pass
        else:
            raise EquationException("n_cold: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        """
        Verify that the prescribed has a valid format.
        """
        if len(self.density.shape) != 2:
            raise EquationException("n_cold: Invalid number of dimensions in prescribed data. Expected 2 dimensions (time x radius).")
        elif len(self.times.shape) != 1:
            raise EquationException("n_cold: Invalid number of dimensions in time grid of prescribed data. Expected one dimension.")
        elif len(self.radius.shape) != 1:
            raise EquationException("n_cold: Invalid number of dimensions in radial grid of prescribed data. Expected one dimension.")
        elif self.density.shape[0] != self.times.size or self.density.shape[1] != self.radius.size:
            raise EquationException("n_cold: Invalid dimensions of prescribed data: {}x{}. Expected {}x{} (time x radius)."
                .format(self.density.shape[0], self.density.shape[1], self.times.size, self.radius.size))


