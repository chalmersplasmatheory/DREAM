# Settings for the runaway electron density

import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity

class RunawayElectrons(UnknownQuantity):
    

    def __init__(self, settings, avalanche=True):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.avalanche = avalanche


    def setAvalanche(self, avalanche):
        """
        Enables/disables avalanche generation.
        """
        self.avalanche = avalanche


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.avalanche = (data['avalanche'] != 0)


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this RunawayElectrons object.
        """
        data = {
            'avalanche': self.avalanche
        }
        
        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if type(self.avalanche) != bool and type(self.avalanche) != int:
            raise EquationException("n_re: Invalid value assigned to 'avalanche'. Expected bool.")


