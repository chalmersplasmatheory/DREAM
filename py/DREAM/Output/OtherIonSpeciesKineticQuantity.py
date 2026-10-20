
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


    def __getitem__(self, name):
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

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None, description=None, momentumgrid=None):
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
        if momentumgrid is None:
            momentumgrid = self.momentumgrid

        return type(self)(name=name, data=data, grid=grid, output=output, momentumgrid=momentumgrid, description=description)

