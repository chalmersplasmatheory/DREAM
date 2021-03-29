# Base class for "other" kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np

from . OutputException import OutputException
from . FluidQuantity import FluidQuantity


class OtherFluidQuantity(FluidQuantity):
    

    def __init__(self, name, data, description, grid, output):
        """
        Constructor.
        """
        attr = {'description': description}
        super(OtherFluidQuantity, self).__init__(name=name, data=data, grid=grid, attr=attr, output=output)

        self.time = grid.t[1:]


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        #s = self.__str__() 
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Other fluid quantity of size NT x NR = {} x {}'.format(self.name, self.data.shape[0], self.data.shape[1])


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def _renormalizeTimeIndexForUnknown(self, t):
        """
        Tries to re-normalize the given time index so that it correctly indexes
        a regular unknown quantity (which has a different time base).
        """
        if t is None:
            t = slice(None)

        start, stop, step = t.start, t.stop, t.step
        if start is None or start == 0:
            start = 1
        if stop is not None:
            stop += 1

        t = slice(start, stop, step)

        return t


