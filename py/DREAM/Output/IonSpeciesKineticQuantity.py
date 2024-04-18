
from . KineticQuantity import KineticQuantity
from . UnknownQuantity import UnknownQuantity


class IonSpeciesKineticQuantity(UnknownQuantity):
    

    def __init__(self, name, data, attr, grid, output, momentumgrid):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)

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


