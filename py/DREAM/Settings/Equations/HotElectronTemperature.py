
import numpy as np
from .EquationException import EquationException
from .ColdElectronTemperature import ColdElectronTemperature
from ..EquationTrigger import EquationTrigger

from .ColdElectronTemperature import TYPE_PRESCRIBED as T_COLD_PRESCRIBED
from .ColdElectronTemperature import HALO_REGION_LOSSES_INCLUDED, HALO_REGION_LOSSES_NEGLECTED
from .ColdElectronTemperature import RECOMBINATION_RADIATION_INCLUDED, RECOMBINATION_RADIATION_NEGLECTED


TYPE_MOMENT = 1
TYPE_SELFCONSISTENT = 2


class HotElectronTemperature(ColdElectronTemperature):
    

    def __init__(
        self, settings, ttype=TYPE_MOMENT,
        temperature=None, radius=0, times=0,
        recombination=RECOMBINATION_RADIATION_NEGLECTED,
        halo_region_losses=HALO_REGION_LOSSES_NEGLECTED,
        makeTrigger=True
    ):
        """
        Constructor.
        """
        super().__init__(
            settings=settings, ttype=T_COLD_PRESCRIBED, # use dummy T_cold value here to avoid problems
            temperature=temperature, radius=radius, times=times,
            recombination=recombination, halo_region_losses=halo_region_losses,
            name='T_hot', makeTrigger=False
        )
        self.setType(ttype)
        self.enabled = False

        if makeTrigger:
            self.trigger = EquationTrigger(
                settings,
                HotElectronTemperature(
                    settings, ttype=ttype, temperature=temperature,
                    radius=radius, times=times, recombination=recombination,
                    halo_region_losses=halo_region_losses, makeTrigger=False
                )
            )
        else:
            self.trigger = None


    def setEnabled(self, enabled=True):
        """
        Enable/disable this quantity in the simulation.
        """
        self.enabled = enabled


    def setType(self, ttype):
        """
        Specifies whether to evolve the hot electron temperature
        self-consistently or as a moment of the distribution function.

        :param ttype: Type of evolution.
        """
        if ttype not in [TYPE_MOMENT, TYPE_SELFCONSISTENT]:
            raise EquationException(f"T_hot: Unrecognized type of equation: {ttype}.")

        self.type = ttype

        if ttype == TYPE_SELFCONSISTENT:
            self.enabled = True


    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        if 'enabled' in data:
            self.enabled = bool(data['enabled'])

        if self.enabled:
            super().fromdict(data)


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this HotElectronTemperature object.
        """
        print(f'T_hot: {self.enabled}')
        if self.enabled:
            data = super().todict()
            data['enabled'] = self.enabled

            if self.type == TYPE_MOMENT:
                del data['data']

            return data
        else:
            return {}


    def verifySettings(self):
        if self.enabled:
            super().verifySettings()


