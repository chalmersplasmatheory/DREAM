# Common class for equation trigger settings

from .. DREAMException import DREAMException


EQN_TRIGGER_TYPE_NONE = 1
EQN_TRIGGER_TYPE_TIME = 2
EQN_TRIGGER_TYPE_COLD_ELECTRON_RISE = 3

ALL_TRIGGER_TYPES = [
    EQN_TRIGGER_TYPE_NONE,
    EQN_TRIGGER_TYPE_TIME,
    EQN_TRIGGER_TYPE_COLD_ELECTRON_RISE
]

class EquationTrigger:
    

    def __init__(self, settings, eqnsettings):
        """
        Constructor.
        """
        self.settings = settings
        self.equation = eqnsettings

        self.sensitivity = 0.01

        # EQN_TRIGGER_TYPE_TIME
        self.trigger_time = 1

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
        if condition not in ALL_TRIGGER_TYPES:
            raise DREAMException(f"Unrecognized trigger condition specified: {condition}.")

        self.condition = condition


    def setTimeTrigger(self, time):
        """
        Enable an equation trigger which switches the equation after a
        given time.
        """
        self.setCondition(EQN_TRIGGER_TYPE_TIME)
        self.trigger_time = time


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
        data = {
            'condition': self.condition,
            'equation': self.equation.todict()
        }

        if self.condition == EQN_TRIGGER_TYPE_TIME:
            data['trigger_time'] = self.trigger_time
        elif self.condition == EQN_TRIGGER_TYPE_COLD_ELECTRON_RISE:
            data['sensitivity'] = self.sensitivity

        return data


