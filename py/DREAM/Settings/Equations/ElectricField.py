
import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . PrescribedInitialParameter import PrescribedInitialParameter
from . PrescribedScalarParameter import PrescribedScalarParameter
from . UnknownQuantity import UnknownQuantity


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2
TYPE_PRESCRIBED_OHMIC_CURRENT = 3

BC_TYPE_PRESCRIBED = 1
BC_TYPE_SELFCONSISTENT = 2
BC_TYPE_TRANSFORMER = 3



class ElectricField(PrescribedParameter, PrescribedInitialParameter, PrescribedScalarParameter, UnknownQuantity):
    
    def __init__(self, settings, ttype=TYPE_PRESCRIBED, efield=None, radius=0,
                 times=0):
        """
        Constructor.

        :param DREAM.DREAMSettings settings: Parent settings object.
        :param int ttype: Method to use for evolving electric field.
        :param efield: Prescribed or initial electric field (profile).
        :param radius: Radial grid on which the prescribed or initial electric field (profile) is defined.
        :param times: Time grid on which the prescribed electric field is defined.
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
        self.bctype = BC_TYPE_PRESCRIBED
        self.inverse_wall_time = None
        self.V_loop_wall_R0 = None
        self.V_loop_wall_t = None
        self.R0 = 0

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
        """
        When ``TYPE_SELFCONSISTENT``, sets the initial electric field profile.
        The parameter ``efield`` may be either a scalar (in which case the
        profile is taken to be uniform) or an array. The associated radial grid
        ``radius`` must be of the same type and dimension.

        :param efield: Initial electric field profile.
        :param radius: Radial grid on which the initial profile is defined.
        """
        _data, _rad = self._setInitialData(data=efield, radius=radius)

        self.efield = _data
        self.radius = _rad
        self.times  = None

        self._verifySettingsPrescribedInitialData()


    def setPrescribedData(self, efield, radius=0, times=0):
        """
        When ``TYPE_PRESCRIBED``, sets the spatiotemporal evolution of the
        electric field during the simulation. The parameter ``efield`` may be
        either a scalar (in which case the electric field is taken to be
        constant and uniform in time and radius) or a 2D array of shape
        (nt, nr). The associated time grid ``times`` must be of size ``nt`` and
        the radial grid must be of size ``nr``.

        :param efield: Prescribed electric field.
        :param radius: Radial grid on which the electric field is prescribed.
        :param times:  Time grid on which the electric field is prescribed.
        """
        _data, _rad, _tim = self._setPrescribedData(efield, radius, times)
        self.efield = _data
        self.radius = _rad
        self.times  = _tim

        self._verifySettingsPrescribedData()


    def setBoundaryCondition(self, bctype = BC_TYPE_SELFCONSISTENT, V_loop_wall_R0=None,
                             times=0, inverse_wall_time=None, R0=0):
        r"""
        Specifies the boundary condition to use when solving for the electric
        field self-consistently, i.e. with ``TYPE_SELFCONSISTENT``. Possible
        boundary condition types are:

        +------------------------+------------------------------------------------------------------------------------------+
        | Name                   | Description                                                                              |
        +========================+==========================================================================================+
        | BC_TYPE_PRESCRIBED     | Set :math:`V_{\rm loop}` on the tokamak wall.                                            |
        +------------------------+------------------------------------------------------------------------------------------+
        | BC_TYPE_SELFCONSISTENT | Specify the tokamak wall time and solve self-consistently for :math:`V_{\rm loop,wall}`. |
        +------------------------+------------------------------------------------------------------------------------------+
        | BC_TYPE_TRANSFORMER    | Same as ``BC_TYPE_SELFCONSISTENT``, but with prescribed loop voltage via transformer.    |
        +------------------------+------------------------------------------------------------------------------------------+

        :param int bctype:        Type of boundary condition to use (see table above for available options).
        :param V_loop_wall_R0:    Prescribed value of :math:`V_{\rm loop}/R_0` on the tokamak wall (or at transformer in case of ``type=BC_TYPE_TRANSFORMER``), normalized to the tokamak major radius :math:`R_0`.
        :param times:             Time grid on which ``V_loop_wall_R0`` is given.
        :param inverse_wall_time: Inverse wall time for the tokamak, used when solving for :math:`V_{\rm loop,wall}` self-consistently.
        :param R0:                Major radius for the tokamak, only used when solving for :math:`V_{\rm loop,wall}` self-consistently (independent of radial-grid major radius).
        """
        if bctype == BC_TYPE_PRESCRIBED:
            self.bctype = bctype

            # Ensure correct format
            _data, _tim = self._setScalarData(data=V_loop_wall_R0, times=times)
            self.V_loop_wall_R0 = _data
            self.V_loop_wall_t = _tim
        elif bctype == BC_TYPE_SELFCONSISTENT:
            self.bctype = bctype
            self.inverse_wall_time = inverse_wall_time
            self.R0 = R0
        elif bctype == BC_TYPE_TRANSFORMER:
            self.bctype = bctype
            self.inverse_wall_time = inverse_wall_time
            self.R0 = R0

            # Ensure correct format
            _data, _tim = self._setScalarData(data=V_loop_wall_R0, times=times)
            self.V_loop_wall_R0 = _data
            self.V_loop_wall_t = _tim
        else:
            raise EquationException("E_field: Unrecognized boundary condition type: {}".format(bctype))


    def setType(self, ttype):
        r"""
        Set the type of equation to use for evolving the electric field. The
        available types are

        +-------------------------------+----------------------------------------------------------------------------------------------+
        | Name                          | Description                                                                                  |
        +===============================+==============================================================================================+
        | TYPE_PRESCRIBED               | Prescribe spatiotemporal evolution of electric field.                                        |
        +-------------------------------+----------------------------------------------------------------------------------------------+
        | TYPE_SELFCONSISTENT           | Evolve electric field consistent with the evolution of the poloidal flux and plasma current. |
        +-------------------------------+----------------------------------------------------------------------------------------------+
        | TYPE_PRESCRIBED_OHMIC_CURRENT | Evolve electric field consistent with the evolution of the poloidal flux and plasma current. |
        +-------------------------------+----------------------------------------------------------------------------------------------+

        :param int ttype: Type of electric field evolution to use.
        """
        if ttype in [TYPE_PRESCRIBED, TYPE_PRESCRIBED_OHMIC_CURRENT]:
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

        if self.type == TYPE_PRESCRIBED:
            self.efield = data['data']['x']
            self.radius = data['data']['r']
            self.times  = data['data']['t']
        elif self.type == TYPE_SELFCONSISTENT:
            self.efield = data['init']['x']
            self.radius = data['init']['r']
            self.bctype = data['bc']['type']
            if self.bctype == BC_TYPE_PRESCRIBED:
                self.V_loop_wall_R0   = data['bc']['V_loop_wall']['x']
                self.V_loop_wall_t = data['bc']['V_loop_wall']['t']
            elif self.bctype == BC_TYPE_SELFCONSISTENT:
                self.inverse_wall_time = data['bc']['inverse_wall_time']
                if not np.isscalar(self.inverse_wall_time):
                    self.inverse_wall_time = float(self.inverse_wall_time[0])
                if 'R0' in data['bc']:
                    self.R0 = data['bc']['R0']
            elif self.bctype == BC_TYPE_TRANSFORMER:
                self.inverse_wall_time = data['bc']['inverse_wall_time']
                if not np.isscalar(self.inverse_wall_time):
                    self.inverse_wall_time = float(self.inverse_wall_time[0])
                if 'R0' in data['bc']:
                    self.R0 = data['bc']['R0']

                self.V_loop_wall_R0 = data['bc']['V_loop_wall']['x']
                self.V_loop_wall_t  = data['bc']['V_loop_wall']['t']
            else:
                raise EquationException("E_field: Unrecognized boundary condition type: {}".format(self.bctype))
        elif self.type == TYPE_PRESCRIBED_OHMIC_CURRENT:
            pass
        else:
            raise EquationException("E_field: Unrecognized electric field type: {}".format(self.type))

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }
        data['bc'] = {'type' : self.bctype}
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
            if self.bctype == BC_TYPE_PRESCRIBED:
                data['bc']['V_loop_wall'] = {
                        'x': self.V_loop_wall_R0,
                        't': self.V_loop_wall_t
                }                
            elif self.bctype == BC_TYPE_SELFCONSISTENT:
                data['bc']['inverse_wall_time'] = self.inverse_wall_time
                data['bc']['R0'] = self.R0
            elif self.bctype == BC_TYPE_TRANSFORMER:
                data['bc']['inverse_wall_time'] = self.inverse_wall_time
                data['bc']['R0'] = self.R0
                data['bc']['V_loop_wall'] = {
                        'x': self.V_loop_wall_R0,
                        't': self.V_loop_wall_t
                }                
        elif self.type == TYPE_PRESCRIBED_OHMIC_CURRENT:
            pass
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

            self._verifySettingsPrescribedData()
        elif self.type == TYPE_SELFCONSISTENT:
            if type(self.efield) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no electric field data provided.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("E_field: Electric field prescribed, but no radial data provided, or provided in an invalid format.")

            # Check boundary condition
            if self.bctype == BC_TYPE_PRESCRIBED:
                self._verifySettingsPrescribedScalarData()
            elif self.bctype == BC_TYPE_SELFCONSISTENT:
                if not np.isscalar(self.inverse_wall_time):
                    raise EquationException("E_field: The specified inverse wall time is not a scalar: {}".format(self.inverse_wall_time))
                if not np.isscalar(self.R0) and not self.R0<0:
                    raise EquationException("E_field: The specified major radius must be scalar and non-negative: R0 = {}".format(self.R0))
            elif self.bctype == BC_TYPE_TRANSFORMER:
                self._verifySettingsPrescribedScalarData()
                if not np.isscalar(self.inverse_wall_time):
                    raise EquationException("E_field: The specified inverse wall time is not a scalar: {}".format(self.inverse_wall_time))
                if not np.isscalar(self.R0) and not self.R0<0:
                    raise EquationException("E_field: The specified major radius must be scalar and non-negative: R0 = {}".format(self.R0))
            else:
                raise EquationException("E_field: Unrecognized boundary condition type: {}.".format(self.bctype))

            self._verifySettingsPrescribedInitialData()
        elif self.type == TYPE_PRESCRIBED_OHMIC_CURRENT:
            pass
        else:
            raise EquationException("E_field: Unrecognized equation type specified: {}.".format(self.type))


    def _verifySettingsPrescribedData(self):
        super()._verifySettingsPrescribedData('E_field', data=self.efield, radius=self.radius, times=self.times)


    def _verifySettingsPrescribedInitialData(self):
        super()._verifySettingsPrescribedInitialData('E_field', data=self.efield, radius=self.radius)


    def _verifySettingsPrescribedScalarData(self):
        super()._verifySettingsPrescribedScalarData('E_field', data=self.V_loop_wall_R0, times=self.V_loop_wall_t)


