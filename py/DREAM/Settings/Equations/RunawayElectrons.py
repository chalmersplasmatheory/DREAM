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

AVALANCHE_MODE_NEGLECT = 1
AVALANCHE_MODE_FLUID = 2
AVALANCHE_MODE_KINETIC = 3

class RunawayElectrons(UnknownQuantity):
    

    def __init__(self, settings, avalanche=AVALANCHE_MODE_NEGLECT, dreicer=DREICER_RATE_DISABLED, Eceff=COLLQTY_ECEFF_MODE_CYLINDRICAL, pCutAvalanche=0):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.avalanche = avalanche
        self.dreicer   = dreicer
        self.Eceff     = Eceff
        self.pCutAvalanche = pCutAvalanche

    def setAvalanche(self, avalanche, pCutAvalanche=0):
        """
        Enables/disables avalanche generation.
        """
        self.avalanche = avalanche
        self.pCutAvalanche = pCutAvalanche


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
        self.avalanche = data['avalanche']
        self.pCutAvalanche = data['pCutAvalanche']
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
            'Eceff': self.Eceff,
            'pCutAvalanche': self.pCutAvalanche
        }
        
        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if type(self.avalanche) != int:
            raise EquationException("n_re: Invalid value assigned to 'avalanche'. Expected integer.")
        if type(self.dreicer) != int:
            raise EquationException("n_re: Invalid value assigned to 'dreicer'. Expected integer.")
        if type(self.Eceff) != int:
            raise EquationException("n_re: Invalid value assigned to 'Eceff'. Expected integer.")
        if self.avalanche == AVALANCHE_MODE_KINETIC and self.pCutAvalanche is 0:
            raise EquationException("n_re: Invalid value assigned to 'pCutAvalanche'. Must be set explicitly when using KINETIC avalanche.")


