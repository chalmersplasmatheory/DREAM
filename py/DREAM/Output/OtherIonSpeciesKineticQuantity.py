
import numpy as np

from. IonsKineticQuantity import IonsKineticQuantity
from . OtherKineticQuantity import OtherKineticQuantity


class OtherIonSpeciesKineticQuantity(IonsKineticQuantity):
    

    def __init__(self, name, data, description, grid, output, momentumgrid):
        """
        Constructor.
        """
        attr = {'description': description}
        super().__init__(name=name, data=data, grid=grid, attr=attr, output=output, momentumgrid=momentumgrid)

        self.time = grid.t[1:]


    def __repr__(self):
        return self.__str__()


    def __getitem(self, name):
        """
        Direct access to data.
        """
        idx = self.ions.getIndex(name)

        return OtherKineticQuantity(
            name=f'{self.name}_{name}',
            data=self.data[:,idx,:],
            description=self.description,
            grid=self.grid,
            output=self.output,
            momentumgrid=self.momentumgrid
        )


