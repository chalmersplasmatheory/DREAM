
import matplotlib.pyplot as plt
import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2


class ColdElectronTemperature(PrescribedParameter):
    
    def __init__(self, ttype=1, temperature=None, radius=None, times=None):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        self.temperature = None
        self.radius  = None
        self.times   = None

        if (ttype == TYPE_PRESCRIBED) and (temperature is not None) and (radius is not None) and (times is not None):
            self.setPrescribedData(temperature, radius=radius, times=times)


    ###################
    # SETTERS
    ###################
    def setPrescribedData(self, temperature, radius=0, times=0):
        _t, _rad, _tim = self._setPrescribedData(temperature, radius, times)
        self.temperature = _t
        self.radius      = _rad
        self.times       = _tim

        self.verifySettingsPrescribedData()


    def setType(self, ttype):
        if ttype == TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == TYPE_SELFCONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))


    def fromdict(self, data):
        self.type = data['type']

        if self.type == TYPE_PRESCRIBED:
            self.temperature = data['data']['x']
            self.radius = data['data']['r']
            self.times = data['data']['t']
        elif self.type == TYPE_SELFCONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }

        if self.type == TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.temperature,
                'r': self.radius,
                't': self.times
            }
        elif self.type == TYPE_SELFCONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.type == TYPE_PRESCRIBED:
            if type(self.temperature) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no temperature data provided.")
            elif type(self.times) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no time data provided, or provided in an invalid format.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedData()
        elif self.type == TYPE_SELFCONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        self._verifySettingsPrescribedData('T_cold', self.temperature, self.radius, self.times)


