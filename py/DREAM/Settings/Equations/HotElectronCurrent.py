
import numpy as np
from .EquationException import EquationException
from .UnknownQuantity import UnknownQuantity
from .. EquationTrigger import EquationTrigger


JHOT_TYPE_MOMENT = 1
JHOT_TYPE_OHMIC = 2


class HotElectronCurrent(UnknownQuantity):
    

    def __init__(self, settings, makeTrigger=True):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.eqn_type = JHOT_TYPE_MOMENT

        if makeTrigger:
            self.trigger = EquationTrigger(
                settings,
                HotElectronCurrent(settings, makeTrigger=False)
            )
        else:
            self.trigger = None


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


