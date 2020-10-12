
import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . PrescribedInitialParameter import PrescribedInitialParameter
from . PrescribedScalarParameter import PrescribedScalarParameter
from . UnknownQuantity import UnknownQuantity


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2

BC_TYPE_PRESCRIBED = 1
BC_TYPE_SELFCONSISTENT = 2



class ElectricField(PrescribedParameter, PrescribedInitialParameter, PrescribedScalarParameter, UnknownQuantity):
    
    def __init__(self, settings, ttype=TYPE_PRESCRIBED, efield=None, radius=0, times=0, wall_radius=-1):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.setType(ttype=ttype)

        # Prescribed electric field evolution
        self.efield = None
        self.radius = None
        self.times  = None

        if (ttype == TYPE_PRESCRIBED) and (efield is not None):
            self.setPrescribedData(efield=efield, radius=radius, times=times)
        elif ttype == TYPE_SELFCONSISTENT:
            self.setInitialProfile(efield=efield, radius=radius)

        # Boundary condition quantities
        self.bctype = None
        self.inverse_wall_time = None
        self.V_loop_wall = None
        self.V_loop_wall_t = None
        self.wall_radius = wall_radius

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


    def setBoundaryCondition(self, bctype = BC_TYPE_SELFCONSISTENT, V_loop_wall=None, times=0, inverse_wall_time=None, wall_radius=-1):
        self.wall_radius = wall_radius
        if bctype == BC_TYPE_PRESCRIBED:
            self.bctype = bctype

            # Ensure correct format
            _data, _tim = self._setScalarData(data=V_loop_wall, times=times)
            self.V_loop_wall = _data
            self.V_loop_wall_t = _tim
        elif bctype == BC_TYPE_SELFCONSISTENT:
            self.bctype = bctype
            self.inverse_wall_time = inverse_wall_time
        else:
            raise EquationException("E_field: Unrecognized boundary condition type: {}".format(bctype))


    def setType(self, ttype):
        if ttype == TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == TYPE_SELFCONSISTENT:
            self.type = ttype

            # Set E=0 if 'setInitialProfile' has not been previously called
            # (if 'setInitialProfile()' has been called, 'self.radius != None'
            # and 'self.times == None')
            if (self.radius) is None or (self.times is not None):
                self.setInitialProfile(efield=0)
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(ttype))


    def fromdict(self, data):
        """
        Sets this paramater from settings provided in a dictionary.
        """
        self.type = data['type']
        self.wall_radius = data['bc']['wall_radius']

        if type(self.wall_radius) == np.ndarray:
            self.wall_radius = float(self.wall_radius[0])
        else:
            self.wall_radius = float(self.wall_radius)

        if self.type == TYPE_PRESCRIBED:
            self.efield = data['data']['x']
            self.radius = data['data']['r']
            self.times  = data['data']['t']
        elif self.type == TYPE_SELFCONSISTENT:
            self.efield = data['init']['x']
            self.radius = data['init']['r']
            self.bctype = data['bc']['type']
            if self.bctype == BC_TYPE_PRESCRIBED:
                self.V_loop_wall   = data['bc']['V_loop_wall']['x']
                self.V_loop_wall_t = data['bc']['V_loop_wall']['t']
            elif self.bctype == BC_TYPE_SELFCONSISTENT:
                self.inverse_wall_time = data['bc']['inverse_wall_time']
                if not np.isscalar(self.inverse_wall_time):
                    self.inverse_wall_time = float(self.inverse_wall_time[0])
            else:
                raise EquationException("E_field: Unrecognized boundary condition type: {}".format(self.bctype))

            
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }
        data['bc'] = {'wall_radius' : self.wall_radius }

        if self.type == TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.efield,
                'r': self.radius,
                't': self.times
            }
        elif self.type == TYPE_SELFCONSISTENT:
            data['init'] = {
                'x': self.efield,
                'r': self.radius,
            }
            data['bc']['type'] = self.bctype
            
            if self.bctype == BC_TYPE_PRESCRIBED:
                data['bc']['V_loop_wall'] = {
                        'x': self.V_loop_wall,
                        't': self.V_loop_wall_t
                }                
            elif self.bctype == BC_TYPE_SELFCONSISTENT:
                data['bc']['inverse_wall_time'] = self.inverse_wall_time
                    
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if not np.isscalar(self.wall_radius):
            raise EquationException("E_field: The specified wall radius is not a scalar: {}.".format(self.wall_radius))
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

            # Check boundary condition
            if self.bctype == BC_TYPE_PRESCRIBED:
                self.verifySettingsPrescribedScalarData()
            elif self.bctype == BC_TYPE_SELFCONSISTENT:
                if not np.isscalar(self.inverse_wall_time):
                    raise EquationException("E_field: The specified inverse wall time is not a scalar: {}".format(self.inverse_wall_time))
            else:
                raise EquationException("E_field: Unrecognized boundary condition type: {}.".format(self.bctype))

            self.verifySettingsPrescribedInitialData()
        else:
            raise EquationException("E_field: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        self._verifySettingsPrescribedData('E_field', data=self.efield, radius=self.radius, times=self.times)


    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('E_field', data=self.efield, radius=self.radius)


    def verifySettingsPrescribedScalarData(self):
        self._verifySettingsPrescribedScalarData('E_field', data=self.V_loop_wall, times=self.V_loop_wall_t)


