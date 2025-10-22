# Common class for equation trigger settings

from .. DREAMException import DREAMException


EQN_TRIGGER_TYPE_NONE = 1
EQN_TRIGGER_TYPE_COLD_ELECTRON_RISE = 2


class EquationTrigger:
    

    def __init__(self, settings, eqnsettings):
        """
        Constructor.
        """
        self.settings = settings
        self.equation = eqnsettings

        self.sensitivity = 0.01

        self.condition = EQN_TRIGGER_TYPE_NONE


    def enabled(self):
        """
        Returns ``True`` if the trigger condition is enabled,
        otherwise ``False``.
        """
        return (self.condition != EQN_TRIGGER_TYPE_NONE)


    def setCondition(self, condition):
        """
        Set the condition to use for triggering an equation switch.
        """
        if condition not in [EQN_TRIGGER_TYPE_NONE, EQN_TRIGGER_TYPE_COLD_ELECTRON_RISE]:
            raise DREAMException(f"Unrecognized trigger condition specified: {condition}.")

        self.condition = condition


    def fromdict(self, data):
        """
        Load settings from a dictionary.
        """
        self.condition = data['condition']
        self.equation.fromdict(data['equation'])


    def todict(self):
        """
        Turn these settings into a dictionary.
        """
        return {
            'condition': self.condition
            'equation': self.equation.todict()
        }


