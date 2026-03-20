

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

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None, description=None):
        """
        Creates a new object of the same type where the provided quantities replace
        those of self.
        """
        if name is None:
            name = self.name
        if data is None:
            data = self.data
        if grid is None:
            grid = self.grid
        if output is None:
            output = self.output
        if description is None:
            if attr is not None and "description" in attr:
                description = attr["description"]
            else:
                description = self.description
        
        return type(self)(name=name, data=data, grid=grid, output=output, description=description)
    
