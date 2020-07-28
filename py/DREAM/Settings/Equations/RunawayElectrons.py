# Settings for the runaway electron density

import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity


DREICER_RATE_DISABLED = 1
DREICER_RATE_CONNOR_HASTIE_NOCORR= 2
DREICER_RATE_CONNOR_HASTIE = 3
DREICER_RATE_NEURAL_NETWORK = 4

COLLQTY_ECEFF_MODE_CYLINDRICAL = 1
COLLQTY_ECEFF_MODE_SIMPLE = 2
COLLQTY_ECEFF_MODE_FULL = 3


class RunawayElectrons(UnknownQuantity):
    

    def __init__(self, settings, avalanche=True, dreicer=DREICER_RATE_DISABLED, Eceff=COLLQTY_ECEFF_MODE_CYLINDRICAL):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.avalanche = avalanche
        self.dreicer   = dreicer
        self.Eceff     = Eceff


    def setAvalanche(self, avalanche):
        """
        Enables/disables avalanche generation.
        """
        self.avalanche = avalanche


    def setDreicer(self, dreicer):
        """
        Specifies which model to use for calculating the
        Dreicer runaway rate.
        """
        self.dreicer = dreicer

    def setEceff(self, Eceff):
        """
        Specifies which model to use for calculating the
        effective critical field (used in the avalanche formula).
        """
        self.Eceff = Eceff


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.avalanche = (data['avalanche'] != 0)
        self.dreicer   = data['dreicer']
        self.Eceff     = data['Eceff']

    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this RunawayElectrons object.
        """
        data = {
            'avalanche': self.avalanche,
            'dreicer': self.dreicer,
            'Eceff': self.Eceff
        }
        
        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if type(self.avalanche) != bool and type(self.avalanche) != int:
            raise EquationException("n_re: Invalid value assigned to 'avalanche'. Expected bool.")
        if type(self.dreicer) != int:
            raise EquationException("n_re: Invalid value assigned to 'dreicer'. Expected integer.")
        if type(self.Eceff) != int:
            raise EquationException("n_re: Invalid value assigned to 'Eceff'. Expected integer.")


