# Representation of a kinetic quantity which has an ion index.

import matplotlib.pyplot as plt
from . IonSpeciesKineticQuantity import IonSpeciesKineticQuantity
from . KineticQuantity import KineticQuantity
from . UnknownQuantity import UnknownQuantity


class IonsKineticQuantity(UnknownQuantity):
    

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
        Convert this object to an "official" string.
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
            nt, ni, nr, nxi, np = self.data.shape
        else:
            ni, nt, nr, nxi, np = 1, *self.data.shape

        s = f'({self.name}) Ion species kinetic quantity of size NI x NT x NR x NXI x NP = {ni} x {nt} x {nr} x {nxi} x {np}\n'
        for i in range(len(self.ions.Z)):
            s += f"  {self.ions[i][0]:2s} (Z = {self.ions[i][1]:3d})\n"

        return s


    def __getitem__(self, name):
        """
        Direct access to data.
        """
        idx = self.ions.getIndex(name)
        offs = self.ions.getIonOffset(name)
        Z = self.ions.getCharge(idx)

        if self.data.ndim == 5:
            data = self.data[:,offs:(offs+Z+1),:]
            return IonSpeciesKineticQuantity(name=f'{self.name}_{name}', data=data, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid, attr=self.attr)
        else:
            data = self.data[:]
            return KineticQuantity(name=f'{self.name}_{name}', data=data, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid, attr=self.attr)


    def dumps(self, ion=None, r=None, t=None, xi=None, p=None):
        """
        Print the data in this quantity.
        """
        return self.get(ion=ion, r=r, t=t, xi=xi, p=p).__str__()


    def get(self, ion=None, r=None, t=None, xi=None, p=None):
        """
        Returns the data for the specified ion, or in the specified time
        interval or radial point. If none of the indices are given, returns
        the full evolution of the quantity.
        """
        sion = ion if ion is not None else slice(None)
        sr   = r if r is not None else slice(None)
        st   = t if t is not None else slice(None)
        sxi  = xi if xi is not None else slice(None)
        sp   = p if p is not None else slice(None)

        if self.data.ndim == 5:
            return self.data[st,sion,sr,sxi,sp]
        else:
            return self.data[st,sr,sxi,sp]


