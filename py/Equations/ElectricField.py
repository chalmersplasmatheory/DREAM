
import numpy as np
from Equations.EquationException import EquationException


class ElectricField:
    
    TYPE_PRESCRIBED = 1

    def __init__(self, ttype=1, efield=None, radius=None, times=None):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Prescribed electric field evolution
        self.efield = None
        self.radius = None
        self.times  = None

        if (ttype == self.TYPE_PRESCRIBED) and (efield is not None) and (radius is not None) and (times is not None):
            self.setPrescribedData(efield=efield, radius=radius, times=times)


    ####################
    # SETTERS
    ####################
    def setPrescribedData(self, efield, radius, times):
        def convtype(v, name):
            if type(v) == list: return np.array(v)
            elif type(v) == np.ndarray: return v
            else: raise EquationException("E_field: Invalid data type of prescribed '{}'.".format(name))

        self.efield  = convtype(efield, 'efield')
        self.radius  = convtype(radius, 'radius')
        self.times   = convtype(times, 'times')

        self.verifySettingsPrescribedData()


    def setType(self, ttype):
        if ttype == self.TYPE_PRESCRIBED:
            self.type = ttype
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }

        if self.type == self.TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.efield,
                'r': self.radius,
                't': self.times
            }
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.type == self.TYPE_PRESCRIBED:
            if type(self.efield) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no electric field data provided.")
            elif type(self.times) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no time data provided, or provided in an invalid format.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedData()
        else:
            raise EquationException("E_field: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        if len(self.efield.shape) != 2:
            raise EquationException("E_field: Invalid number of dimensions in prescribed data. Expected 2 dimensions (time x radius).")
        elif len(self.times.shape) != 1:
            raise EquationException("E_field: Invalid number of dimensions in time grid of prescribed data. Expected one dimension.")
        elif len(self.radius.shape) != 1:
            raise EquationException("E_field: Invalid number of dimensions in radial grid of prescribed data. Expected one dimension.")
        elif self.efield.shape[0] != self.times.size or self.efield.shape[1] != self.radius.size:
            raise EquationException("E_field: Invalid dimensions of prescribed data: {}x{}. Expected {}x{} (time x radius)."
                .format(self.efield.shape[0], self.efield.shape[1], self.times.size, self.radius.size))


