
import matplotlib.pyplot as plt
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException


class ColdElectronTemperature:
    
    TYPE_PRESCRIBED = 1
    TYPE_SELF_CONSISTENT = 2

    def __init__(self, ttype=1, temperature=None, radius=None, times=None):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        self.temperature = None
        self.radius  = None
        self.times   = None

        if (ttype == self.TYPE_PRESCRIBED) and (temperature is not None) and (radius is not None) and (times is not None):
            self.setPrescribedData(temperature, radius=radius, times=times)


    ###################
    # SETTERS
    ###################
    def setPrescribedData(self, temperature, radius, times):
        def convtype(v, name):
            if type(v) == list: return np.array(v)
            elif type(v) == np.ndarray: return v
            else: raise EquationException("T_cold: Invalid data type of prescribed '{}'.".format(name))

        self.temperature = convtype(temperature, 'temperature')
        self.radius  = convtype(radius, 'radius')
        self.times   = convtype(times, 'times')

        self.verifySettingsPrescribedData()


    def setType(self, ttype):
        if ttype == self.TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == self.TYPE_SELF_CONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }

        if self.type == self.TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.temperature,
                'r': self.radius,
                't': self.times
            }
        elif self.type == self.TYPE_SELF_CONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.type == self.TYPE_PRESCRIBED:
            if type(self.temperature) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no temperature data provided.")
            elif type(self.times) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no time data provided, or provided in an invalid format.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedData()
        elif self.type == self.TYPE_SELF_CONSISTENT:
            raise EquationException("T_cold: Self-consistent temperature evolution is not yet supported.")
        else:
            raise EquationException("T_cold: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        if len(self.temperature.shape) != 2:
            raise EquationException("T_cold: Invalid number of dimensions in prescribed data. Expected 2 dimensions (time x radius).")
        elif len(self.times.shape) != 1:
            raise EquationException("T_cold: Invalid number of dimensions in time grid of prescribed data. Expected one dimension.")
        elif len(self.radius.shape) != 1:
            raise EquationException("T_cold: Invalid number of dimensions in radial grid of prescribed data. Expected one dimension.")
        elif self.temperature.shape[0] != self.times.size or self.temperature.shape[1] != self.radius.size:
            raise EquationException("T_cold: Invalid dimensions of prescribed data: {}x{}. Expected {}x{} (time x radius)."
                .format(self.temperature.shape[0], self.temperature.shape[1], self.times.size, self.radius.size))


