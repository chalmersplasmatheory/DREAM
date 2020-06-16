
import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . PrescribedInitialParameter import PrescribedInitialParameter


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2

BC_TYPE_PRESCRIBED = 1
BC_TYPE_SELFCONSISTENT = 2



class ElectricField(PrescribedParameter, PrescribedInitialParameter):
    
    def __init__(self, ttype=TYPE_PRESCRIBED, efield=None, radius=0, times=0, wall_radius=-1):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Prescribed electric field evolution
        self.efield = None
        self.radius = None
        self.times  = None
        self.wall_radius = wall_radius

        
        if (ttype == TYPE_PRESCRIBED) and (efield is not None):
            self.setPrescribedData(efield=efield, radius=radius, times=times)
        elif ttype == TYPE_SELFCONSISTENT:
            self.setInitialProfile(efield=efield, radius=radius)

        # Boundary condition quantities
        self.bctype = None
        self.inverse_wall_time = None
        self.V_loop_wall = None
        self.V_loop_wall_t = None

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

    def setBoundaryCondition(self, bctype = BC_TYPE_SELFCONSISTENT, V_loop_wall=None, times=0, inverse_wall_time=None):
        if bctype == BC_TYPE_PRESCRIBED:
            self.bctype = bctype
            # TODO
        elif bctype == BC_TYPE_SELFCONSISTENT:
            self.bctype = bctype
            self.inverse_wall_time = inverse_wall_time
            self.V_loop_wall = V_loop_wall
            self.V_loop_wall_t = times
        else:
            raise EquationException("E_field: Unrecognized boundary condition type: {}".format(bctype))



    def setType(self, ttype):
        if ttype == TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == TYPE_SELFCONSISTENT:
            self.type = ttype
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(ttype))


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
            self.wall_radius = data['bc']['wall_radius']
            self.inverse_wall_time = data['bc']['inverse_wall_time']
            self.bctype = data['bc']['type']
            self.V_loop_wall   = data['bc']['data']['x']
            self.V_loop_wall_t = data['bc']['data']['t']
            
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
                'r': self.radius,
            }
            # TODO: somehow only write data if bctype is Prescribed
            data['bc'] = {
                'inverse_wall_time': self.inverse_wall_time,
                'wall_radius': self.wall_radius,
                'type': self.bctype,
                'data': {
                    'x': self.V_loop_wall,
                    't': self.V_loop_wall_t
                }
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

        # TODO: verify boundary condition settings 

    def verifySettingsPrescribedData(self):
        self._verifySettingsPrescribedData('E_field', data=self.efield, radius=self.radius, times=self.times)


    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('E_field', data=self.efield, radius=self.radius)


