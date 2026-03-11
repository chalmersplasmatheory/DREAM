

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

        # Handle both 2D (single species) and 3D (multiple species) data
        if self.data.ndim == 3:
            species_data = self.data[:,idx,:]
        elif self.data.ndim == 2:
            # Single species case - data is (NT, NR)
            if idx != 0:
                raise ValueError(f"Species index {idx} out of range for single-species data")
            species_data = self.data
        else:
            raise ValueError(f"Unexpected data dimensions: {self.data.ndim}")

        return OtherFluidQuantity(name='{}_{}'.format(self.name, name), data=species_data, description=self.description, grid=self.grid, output=self.output)

    
