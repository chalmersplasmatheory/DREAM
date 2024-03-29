# Base class for "other" kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np

from . OutputException import OutputException
from . ScalarQuantity import ScalarQuantity


class OtherScalarQuantity(ScalarQuantity):
    

    def __init__(self, name, data, description, grid, output):
        """
        Constructor.
        """
        attr = {'description': description}
        super(OtherScalarQuantity, self).__init__(name=name, data=data, grid=grid, attr=attr, output=output)

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
        return '({}) Other scalar quantity of size NT = {}'.format(self.name, self.data.shape[0])


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

        
