# Representation of a scalar quantity which has an ion index.

from . ScalarQuantity import ScalarQuantity
from . UnknownQuantity import UnknownQuantity


class IonSpeciesScalarQuantity(UnknownQuantity):
    

    def __init__(self, name, data, attr, grid, output):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)

        self.attr = attr
        self.ions = output.ionmeta


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        s = self.__str__()

        if hasattr(self, 'description') and hasattr(self, 'description_eqn'):
            s += "\n:: {}\n:: Evolved using: {}\n".format(self.description, self.description_eqn)

        return s


    def __str__(self):
        """
        Convert this object to a string.
        """
        s = '({}) Ion species scalar quantity of size NI x NT = {} x {}\n'.format(self.name, *self.data.shape)
        for i in range(len(self.ions.Z)):
            s += "  {:2s} (Z = {:3d})\n".format(*self.ions[i])

        return s


    def __getitem__(self, name):
        """
        Direct access to data.
        """
        idx = self.ions.getIndex(name)

        return ScalarQuantity(name='{}_{}'.format(self.name, name), data=self.data[:,idx,:], grid=self.grid, output=self.output, attr=self.attr)

    
    def dumps(self, ion=None, r=None, t=None):
        """
        Print the data in this quantity.
        """
        return self.get(ion=ion, r=r, t=t).__str__()


    def get(self, ion=None, t=None):
        """
        Returns data for the specified ion, or in the specified time
        interval or point. If none of the indices are given, returns
        the full evolution of the quantity.
        """
        sion = ion if ion is not None else slice(None)
        st = t if t is not None else slice(None)

        return self.data[st,sion]


