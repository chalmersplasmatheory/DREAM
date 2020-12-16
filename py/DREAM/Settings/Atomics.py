# Settings related to atomic data  

import numpy as np
from DREAM.DREAMException import DREAMException

# ADAS interpolation method
ADAS_INTERP_BILINEAR = 1
ADAS_INTERP_BICUBIC = 2


class Atomics:
    def __init__(self,
            adas_interpolation = ADAS_INTERP_BICUBIC):
        """
        Constructor.
        """
        self.adas_interpolation = adas_interpolation
    
    def fromdict(self, data):
        """
        Load settings from dictionary.
        """
        self.adas_interpolation = data['adas_interpolation']

    def todict(self, verify=True):
        """
        Returns these settings as a dictionary.
        """
        if verify:
            self.verifySettings()
        data = { 'adas_interpolation' : self.adas_interpolation }
        
        return data
        
    
    def verifySettings(self):
        """
        TODO
        """
        pass
