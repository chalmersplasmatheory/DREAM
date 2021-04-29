# Current density sub-class

import numpy as np

from . FluidQuantity import FluidQuantity


class CurrentDensity(FluidQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super(CurrentDensity, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output)


    def current(self, t=None):
        """
        Calculates the total current corresponding to this current density.
        """
        geom = self.grid.GR0/self.grid.Bmin * self.grid.FSA_R02OverR2

        return self.integral(t=t, w=geom) / (2*np.pi)


    def totalCurrent(self, t=None):
        """
        Calculates the total current corresponding to this current density
        (alias for ``current()``).
        """
        return self.current(t=t)
