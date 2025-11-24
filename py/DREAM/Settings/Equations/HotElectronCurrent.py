
import numpy as np
from .EquationException import EquationException
from .UnknownQuantity import UnknownQuantity
from .. EquationTrigger import EquationTrigger


TYPE_MOMENT = 1
TYPE_OHMIC = 2


class HotElectronCurrent(UnknownQuantity):
    

    def __init__(self, settings, makeTrigger=True):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.eqn_type = TYPE_MOMENT

        if makeTrigger:
            self.trigger = EquationTrigger(
                settings,
                HotElectronCurrent(settings, makeTrigger=False),
                self
            )
        else:
            self.trigger = None


    def setType(self, eqn_type):
        """
        Set the type of equation to use for evolving j_hot.
        """
        if eqn_type in [TYPE_MOMENT, TYPE_OHMIC]:
            self.eqn_type = eqn_type
        else:
            raise EquationException(f"Unrecognized equation type for 'j_hot': {eqn_type}.")


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.eqn_type = int(data['type'])

        if 'switch' in data and self.trigger is not None:
            self.trigger.fromdict(data['switch'])
    

    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this class.
        """
        data = {
            'type': self.eqn_type
        }

        if self.trigger is not None and self.trigger.enabled():
            data['switch'] = self.trigger.todict()

        return data


    def verifySettings(self):
        pass


