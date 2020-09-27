# Wrapper class for 'f_hot' -- the hot electron distribution function

import numpy as np
import matplotlib.pyplot as plt

from . DistributionFunction import DistributionFunction


class HotElectronDistributionFunction(DistributionFunction):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output, momentumgrid=grid.hottail)


    def getTeXName(self):
        """
        Returns the TeX-compatible name of this quantity.
        """
        return r'$f_{\rm hot}$'


