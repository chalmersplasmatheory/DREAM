# Common class for equation trigger settings

from .. DREAMException import DREAMException


TYPE_NONE = 1
TYPE_TIME = 2
TYPE_COLD_ELECTRON_RISE = 3

ALL_TRIGGER_TYPES = [
    TYPE_NONE,
    TYPE_TIME,
    TYPE_COLD_ELECTRON_RISE
]

class EquationTrigger:
    

    def __init__(self, settings, eqnsettings):
        """
        Constructor.
        """
        self.settings = settings
        self.equation = eqnsettings

        self.sensitivity = 0.01

        # TYPE_TIME
        self.trigger_time = 1

        self.condition = TYPE_NONE


    def enabled(self):
        """
        Returns ``True`` if the trigger condition is enabled,
        otherwise ``False``.
        """
        return (self.condition != TYPE_NONE)


    def setCondition(self, condition, trigger_time=None):
        """
        Set the condition to use for triggering an equation switch.
        """
        if condition not in ALL_TRIGGER_TYPES:
            raise DREAMException(f"Unrecognized trigger condition specified: {condition}.")

        self.condition = condition

        if condition == TYPE_TIME:
            self.trigger_time = trigger_time
        elif condition == TYPE_COLD_ELECTRON_RISE:
            pass


    def setTimeTrigger(self, time):
        """
        Enable an equation trigger which switches the equation after a
        given time.
        """
        self.setCondition(TYPE_TIME, trigger_time=time)


    def fromdict(self, data):
        """
        Load settings from a dictionary.
        """
        if 'condition' in data:
            self.condition = data['condition']
        if 'equation' in data:
            self.equation.fromdict(data['equation'])


    def todict(self):
        """
        Turn these settings into a dictionary.
        """
        data = {
            'condition': self.condition,
            'equation': self.equation.todict()
        }

        if self.condition == TYPE_TIME:
            data['trigger_time'] = self.trigger_time
        elif self.condition == TYPE_COLD_ELECTRON_RISE:
            data['sensitivity'] = self.sensitivity

        return data


