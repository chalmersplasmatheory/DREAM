import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings

CONDUCTIVITY_MODE_BRAAMS = 1
CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS = 2
CONDUCTIVITY_MODE_SAUTER_COLLISIONAL = 3

CORRECTED_CONDUCTIVITY_DISABLED = 1
CORRECTED_CONDUCTIVITY_ENABLED  = 2

class OhmicCurrent(UnknownQuantity):
    
    def __init__(self, settings, condMode=CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS, corrCond=CORRECTED_CONDUCTIVITY_ENABLED):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.condMode=condMode
        self.corrCond=corrCond

    def setCorrectedConductivity(self, mode):
        r"""
        Specifies whether to use a conductivity correction, which is added
        to the cold current carried by the distribution function
        
        :param int mode:    Type of conductivity correction to use
        """
        if type(mode) == bool:
            self.corrCond = CORRECTED_CONDUCTIVITY_ENABLED if mode else CORRECTED_CONDUCTIVITY_DISABLED
        else:
            self.corrCond = int(mode) 

    def setConductivityMode(self, mode):
        r"""
        Specifies the formula to use for the condictivity in the ohmic current. 
        The Sauter models are based on O Sauter, C Angioni, YR Lin-Liu, PoP (1999).
        The underlying base 'Spitzer' (slab) conductivity used is the relativistic
        BJ Braams, CFF Karney PoF (1989) formula, where we interpolate in the 
        tabulated values from the paper.

        Possible modes are:

        +----------------------------------------+------------------------------------------------------------------------------------------+
        | Name                                   | Description                                                                              |
        +========================================+==========================================================================================+
        | CONDUCTIVITY_MODE_BRAAMS               | Uses the base Braams & Karney formula                                                    |
        +----------------------------------------+------------------------------------------------------------------------------------------+
        | CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS | Using the Sauter neoclassical correction in the collisionless limit (banana regime).     |
        +----------------------------------------+------------------------------------------------------------------------------------------+
        | CONDUCTIVITY_MODE_SAUTER_COLLISIONAL   | Uses the full Sauter model with collisional neoclassical corrections                     |
        +----------------------------------------+------------------------------------------------------------------------------------------+

        :param int mode:    Type of model to use for the plasma conductivity
        """
        self.condMode = int(mode)


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        if 'conductivityMode' in data:
            self.condMode = data['conductivityMode']
        if 'correctedConductivity' in data:
            self.corrCond = data['correctedConductivity']


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this PoloidalFlux object.
        """
        data = {
            'conductivityMode': self.condMode,
            'correctedConductivity': self.corrCond
        }

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        pass #TODO


