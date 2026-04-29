
from . KineticQuantity import KineticQuantity
from . UnknownQuantity import UnknownQuantity


class IonSpeciesKineticQuantity(UnknownQuantity):
    

    def __init__(self, name, data, attr, grid, output, momentumgrid, triggerinfo=None):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output, triggerinfo=triggerinfo)

        self.attr = attr
        self.ions = output.ionmeta
        self.momentumgrid = momentumgrid


    def __repr__(self):
        """
        Convert this object into an "official" string.
        """
        s = self.__str__()

        if hasattr(self, 'description') and hasattr(self, 'description_eqn'):
            s += f"\n:: {self.description}\n:: Evolved using: {self.description_eqn}\n"

        return s


    def __str__(self):
        """
        Convert this object into a string.
        """
        if self.data.ndim == 5:
            nt, nZ0, nr, nxi, np = self.data.shape
        else:
            nZ0, nt, nr, nxi, np = 1, *self.data.shape

        s = f'({self.name}) Ion species kinetic quantity of size NZ0 x NT x NR x NXI x NP = {nZ0} x {nt} x {nr} x {nxi} x {np}\n'

        return s


    def __getitem__(self, Z0):
        """
        Direct access to data.
        """
        if self.data.ndim == 5:
            data = self.data[:,Z0,:]
        else:
            data = self.data[:]

        return KineticQuantity(name=f'{self.name}_{Z0}', data=data, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid, attr=self.attr)

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None, momentumgrid=None):
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
        if attr is None:
            attr = {}
            if hasattr(self, "description"):
                attr["description"] = self.description
            if hasattr(self, "description_eqn"):
                attr["equation"] = self.description_eqn
        if momentumgrid is None:
            momentumgrid = self.momentumgrid

        return type(self)(name=name, data=data, grid=grid, output=output, momentumgrid=momentumgrid, attr=attr)

