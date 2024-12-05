
import matplotlib.pyplot as plt
import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . PrescribedInitialParameter import PrescribedInitialParameter
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2

RECOMBINATION_RADIATION_INCLUDED = True
RECOMBINATION_RADIATION_NEGLECTED = False

HALO_REGION_LOSSES_INCLUDED = True
HALO_REGION_LOSSES_NEGLECTED = False


class ColdElectronTemperature(PrescribedParameter,PrescribedInitialParameter,UnknownQuantity):
    
    def __init__(self, settings, ttype=TYPE_PRESCRIBED, temperature=None, radius=0, times=0, 
                 recombination=RECOMBINATION_RADIATION_NEGLECTED, halo_region_losses = HALO_REGION_LOSSES_INCLUDED):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.setType(ttype=ttype)

        self.temperature = None
        self.radius = None
        self.times  = None

        self.transport = TransportSettings(kinetic=False)
        self.recombination = recombination

        self.halo_region_losses = halo_region_losses

        if (ttype == TYPE_PRESCRIBED) and (temperature is not None):
            self.setPrescribedData(temperature=temperature, radius=radius, times=times)
        elif ttype == TYPE_SELFCONSISTENT:
            self.setInitialProfile(temperature=temperature, radius=radius)


    ###################
    # SETTERS
    ###################
    def setInitialProfile(self, temperature, radius=0):
        """
        Sets the initial temperature profile T=T(r) for when the temperature is
        evolved self-consistently.

        :param temperature: Scalar or vector giving the initial temperature profile.
        :param radius: If ``temperature`` is a vector, contains the corresponding radial grid on which ``temperature`` is defined.
        """
        _data, _rad = self._setInitialData(data=temperature, radius=radius)

        self.temperature = _data
        self.radius      = _rad
        self.times       = None

        self.verifySettingsPrescribedInitialData()


    def setPrescribedData(self, temperature, radius=0, times=0):
        """
        Prescribes a temperature evolution in time and space.

        :param temperature: Scalar, vector or matrix giving the temperature throughout the simulation.
        :param radius: If ``temperature`` is a function of radius, contains the radial grid on which it is defined.
        :param times: If ``temperature`` is a function of time, contains the time grid on which it is defined.
        """
        _t, _rad, _tim = self._setPrescribedData(temperature, radius, times)
        self.temperature = _t
        self.radius      = _rad
        self.times       = _tim

        self.verifySettingsPrescribedData()


    def setType(self, ttype):
        """
        Specifies whether to evolve the electron temperature according to a
        prescribed function, or self-consistently.

        :param ttype: Type of evolution. Can take one of the following values:

        - ``TYPE_PRESCRIBED``: Evolve according to prescribed function.
        - ``TYPE_SELFCONSISTENT``: Evolve self-consistently.
        """
        if ttype == TYPE_PRESCRIBED:
            self.type = ttype
        elif ttype == TYPE_SELFCONSISTENT:
            self.type = ttype

            # Set T=0 if 'setInitialProfile' has not been previously called
            # (if 'setInitialProfile()' has been called, 'self.radius != None'
            # and 'self.times == None')
            if (self.radius) is None or (self.times is not None):
                self.setInitialProfile(temperature=-1)
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))


    def setRecombinationRadiation(self, recombination=RECOMBINATION_RADIATION_NEGLECTED):
        """
        Specify whether or not to include recombination radiation when evolving
        the temperature self-consistently.
        """
        self.recombination = recombination


    def setHaloRegionLosses(self, halo_region_losses=HALO_REGION_LOSSES_INCLUDED):
        """
        Specify whether or not to include heat losses in the halo region when
        evolving the temperature self-consistently.
        """
        self.halo_region_losses = halo_region_losses


    def fromdict(self, data):
        self.type = data['type']

        if self.type == TYPE_PRESCRIBED:
            self.temperature = data['data']['x']
            self.radius = data['data']['r']
            self.times = data['data']['t']
        elif self.type == TYPE_SELFCONSISTENT:
            self.temperature = data['init']['x']
            self.radius = data['init']['r']

            if 'transport' in data:
                self.transport.fromdict(data['transport'])
            
            if 'halo_region_losses' in data:
                self.halo_region_losses = int(data['halo_region_losses'])
        else:
            raise EquationException("T_cold: Unrecognized cold electron temperature type: {}".format(self.type))
        if 'recombination' in data:
            self.recombination = data['recombination']

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this ColdElectrons object.
        """
        data = { 'type': self.type }
        data['recombination'] = self.recombination
        if self.type == TYPE_PRESCRIBED:
            data['data'] = {
                'x': self.temperature,
                'r': self.radius,
                't': self.times
            }
        elif self.type == TYPE_SELFCONSISTENT:
            data['halo_region_losses'] = self.halo_region_losses
            data['init'] = {
                'x': self.temperature,
                'r': self.radius
            }
            data['transport'] = self.transport.todict()
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
            if type(self.temperature) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no temperature data provided.")
            elif type(self.radius) != np.ndarray:
                raise EquationException("T_cold: Temperature prescribed, but no radial data provided, or provided in an invalid format.")

            self.verifySettingsPrescribedInitialData()
            self.transport.verifySettings()
        else:
            raise EquationException("T_cold: Unrecognized equation type specified: {}.".format(self.type))


    def verifySettingsPrescribedData(self):
        self._verifySettingsPrescribedData('T_cold', self.temperature, self.radius, self.times)

    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('T_cold', data=self.temperature, radius=self.radius)

