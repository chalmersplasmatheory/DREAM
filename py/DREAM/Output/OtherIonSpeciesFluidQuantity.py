

import numpy as np

from . IonSpeciesFluidQuantity import IonSpeciesFluidQuantity
from . OtherFluidQuantity import OtherFluidQuantity


class OtherIonSpeciesFluidQuantity(IonSpeciesFluidQuantity):
    

    def __init__(self, name, data, description, grid, output):
        """
        Constructor.
        """
        attr = {'description': description}
        super().__init__(name=name, data=data, grid=grid, attr=attr, output=output)

        self.time = grid.t[1:]

    
    def __repr__(self):
        return self.__str__()


    def __getitem__(self, name):
        """
        Direct access to data.
        """
        idx = self.ions.getIndex(name)

        return OtherFluidQuantity(name='{}_{}'.format(self.name, name), data=self.data[:,idx,:], description=self.description, grid=self.grid, output=self.output)

    
