# Collision handler

import numpy as np
from DREAM.DREAMException import DREAMException

# Collisionality settings
BREMSSTRAHLUNG_MODE_NEGLECT = 1
BREMSSTRAHLUNG_MODE_STOPPING_POWER = 2
#BREMSSTRAHLUNG_MODE_BOLTZMANN = 3

COLLFREQ_MODE_SUPERTHERMAL = 1
COLLFREQ_MODE_FULL = 2
COLLFREQ_MODE_ULTRA_RELATIVISTIC = 3

COLLFREQ_TYPE_COMPLETELY_SCREENED = 1
COLLFREQ_TYPE_NON_SCREENED = 2
COLLFREQ_TYPE_PARTIALLY_SCREENED = 3

#NONLINEAR_MODE_NEGLECT = 1
#NONLINEAR_MODE_NON_REL_ISOTROPIC = 2
#NONLINEAR_MODE_FULL = 3

LNLAMBDA_CONSTANT = 1
LNLAMBDA_ENERGY_DEPENDENT = 2
LNLAMBDA_THERMAL = 3


PSTAR_MODE_COLLISIONAL = 1
PSTAR_MODE_COLLISIONLESS = 2


class CollisionHandler:
    

    def __init__(self,
            bremsstrahlung_mode=BREMSSTRAHLUNG_MODE_NEGLECT,
            collfreq_mode=COLLFREQ_MODE_FULL, collfreq_type=COLLFREQ_TYPE_PARTIALLY_SCREENED,
            lnlambda=LNLAMBDA_THERMAL, pstar_mode=PSTAR_MODE_COLLISIONLESS):
        """
        Constructor.
        """
        self.bremsstrahlung_mode = bremsstrahlung_mode
        self.collfreq_mode = collfreq_mode
        self.collfreq_type = collfreq_type
        self.lnlambda = lnlambda
        self.pstar_mode = pstar_mode

    
    def fromdict(self, data):
        """
        Load settings from dictionary.
        """
        self.bremsstrahlung_mode = data['bremsstrahlung_mode']
        self.collfreq_mode = data['collfreq_mode']
        self.collfreq_type = data['collfreq_type']
        self.lnlambda = data['lnlambda']
        self.pstar_mode = data['pstar_mode']


    def todict(self, verify=True):
        """
        Returns these settings as a dictionary.
        """
        if verify:
            self.verifySettings()

        return {
            'bremsstrahlung_mode': self.bremsstrahlung_mode,
            'collfreq_mode': self.collfreq_mode,
            'collfreq_type': self.collfreq_type,
            'lnlambda': self.lnlambda,
            'pstar_mode': self.pstar_mode
        }

    
    def verifySettings(self):
        """
        TODO
        """
        pass

