import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings

class PoloidalFlux(UnknownQuantity):
    
    def __init__(self, settings):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.transport = TransportSettings(kinetic=False)

    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        if 'transport' in data:
            self.transport.fromdict(data['transport'])


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this PoloidalFlux object.
        """
        data = {
            'transport': self.transport.todict()
        }

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        self.transport.verifySettings()


