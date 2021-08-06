

import numpy as np

from . IonSpeciesScalarQuantity import IonSpeciesScalarQuantity
from . OtherScalarQuantity import OtherScalarQuantity


class OtherIonSpeciesScalarQuantity(IonSpeciesScalarQuantity):
    

    def __init__(self, name, data, description, grid, output, momentumgrid=None):
        """
        Constructor.
        """
        attr = {'description': description}
        super().__init__(name=name, data=data, grid=grid, attr=attr, output=output)

        self.time = grid.t[1:]

    
    def __repr__(self):
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        s = '({}) Other ion species scalar quantity of size NI x NT = {} x {}\n'.format(self.name, *self.data.shape)
        for i in range(len(self.ions.Z)):
            s += "  {:2s} (Z = {:3d})\n".format(*self.ions[i])

        return s


    def __getitem__(self, name):
        """
        Direct access to data.
        """
        idx = self.ions.getIndex(name)

        return OtherScalarQuantity(name='{}_{}'.format(self.name, name), data=self.data[:,idx,:], description=self.description, grid=self.grid, output=self.output)

    
