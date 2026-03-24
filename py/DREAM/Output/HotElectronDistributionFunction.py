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

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None):
        """
        Creates a new object of the same type where the provided quantities replace
        those of self.

        Unlike KineticQuantity from which this class inherits from,
        the constructor wants no momentumgrid. Therefore overriding that behavior here.
        """
        if name is None:
            name = self.name
        if data is None:
            data = self.data
        if grid is None:
            grid = self.grid
        if output is None:
            output = self.output
        if attr is None:
            attr = {}
            if hasattr(self, "description"):
                attr["description"] = self.description
            if hasattr(self, "description_eqn"):
                attr["equation"] = self.description_eqn

        return type(self)(name=name, data=data, grid=grid, output=output, attr=attr)
