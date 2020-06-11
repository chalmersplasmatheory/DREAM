
import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . PrescribedInitialParameter import PrescribedInitialParameter


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2


class ElectricField(PrescribedParameter, PrescribedInitialParameter):
    
    def __init__(self, ttype=TYPE_PRESCRIBED, efield=None, radius=0, times=0):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Prescribed electric field evolution
        self.efield = None
        self.radius = None
        self.times  = None

        if (ttype == TYPE_PRESCRIBED) and (efield is not None):
            self.setPrescribedData(efield=efield, radius=radius, times=times)
        elif ttype == TYPE_SELFCONSISTENT:
            self.setInitialProfile(efield=efield, radius=radius)


    def __getitem__(self, index):
        """
        Returns the value of the prescribed electric field at
        the given indices.
        """
        return self.efield[index]


    ####################
    # SETTERS
    ####################
    def setInitialProfile(self, efield, radius=0):
        _data, _rad = self._setInitialData(data=efield, radius=radius)

        self.efield = _data
        self.radius = _rad
        self.times  = None

        self.verifySettingsPrescribedInitialData()


    def setPrescribedData(self, efield, radius=0, times=0):
        _data, _rad, _tim = self._setPrescribedData(efield, radius, times)
        self.efield = _data
        self.radius = _rad
        self.times  = _tim

        self.verifySettingsPrescribedData()


    def setType(self, ttype):
        if ttype == TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == TYPE_SELFCONSISTENT:
            self.type = ttype
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))


    def fromdict(self, data):
        """
        Sets this paramater from settings provided in a dictionary.
        """
        self.type = data['type']

        if self.type == TYPE_PRESCRIBED:
            self.efield = data['data']['x']
            self.radius = data['data']['r']
            self.times  = data['data']['t']
        elif self.type == TYPE_SELFCONSISTENT:
            self.efield = data['init']['x']
            self.radius = data['init']['r']
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }

        if self.type == TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.efield,
                'r': self.radius,
                't': self.times
            }
        elif self.type == TYPE_SELFCONSISTENT:
            data['init'] = {
                'x': self.efield,
                'r': self.radius
            }
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.type == TYPE_PRESCRIBED:
            if type(self.efield) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no electric field data provided.")
            elif type(self.times) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no time data provided, or provided in an invalid format.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedData()
        elif self.type == TYPE_SELFCONSISTENT:
            if type(self.efield) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no electric field data provided.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedInitialData()
        else:
            raise EquationException("E_field: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        self._verifySettingsPrescribedData('E_field', data=self.efield, radius=self.radius, times=self.times)


    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('E_field', data=self.efield, radius=self.radius)


