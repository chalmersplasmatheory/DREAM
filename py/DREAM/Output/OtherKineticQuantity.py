# Base class for "other" kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np

from . OutputException import OutputException
from . KineticQuantity import KineticQuantity


class OtherKineticQuantity(KineticQuantity):
    

    def __init__(self, name, data, description, grid, output, momentumgrid=None):
        """
        Constructor.
        """
        attr = {'description': description}
        super(OtherKineticQuantity, self).__init__(name=name, data=data, grid=grid, attr=attr, output=output, momentumgrid=momentumgrid)

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
        return '({}) Other kinetic quantity of size NT x NR x NP2 x NP1 = {} x {} x {} x {}'.format(self.name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3])


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of ion species and
        charge states) covered by this quantity. The total number of elements
        in 'self.data' is the size of the grid on which this quantity lives
        (i.e. scalar grid, fluid grid, or a kinetic grid) times this number.
        """
        return 1

        
