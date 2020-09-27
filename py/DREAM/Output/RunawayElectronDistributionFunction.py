# Wrapper class for 'f_re' -- the runaway electron distribution function

import numpy as np
import matplotlib.pyplot as plt

from .DistributionFunction import DistributionFunction

class RunawayElectronDistributionFunction(DistributionFunction):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output, momentumgrid=grid.runaway)


    def getTeXName(self):
        """
        Returns the TeX-compatible name of this quantity.
        """
        return r'$f_{\rm RE}$'
