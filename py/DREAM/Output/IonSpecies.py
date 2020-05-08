# Object representing a single ion species

import numpy as np
from .IonState import IonState
from .OutputException import OutputException


class IonSpecies:
    

    def __init__(self, name, Z, data, grid, output):
        """
        Constructor.

        name: Name of ion species.
        Z:    Ion atomic charge.
        data: Ion density data (size nt x nZ0 x nr, where nZ0 = Z+1
              is the number of charge states for this ion).
        """
        self.name = name
        self.Z    = Z
        self.grid = grid
        self.output = output
        self.ionstates = list()

        nt = len(self.grid.t)
        nr = len(self.grid.r)

        if data.shape != (nt, Z+1, nr):
            raise OutputException("Invalid dimensions of data. Expected {}x{}x{}, but found {}x{}x{}."
                .format(nt, Z+1, nr, data.shape[0], data.shape[1], data.shape[2]))

        for Z0 in range(0, Z+1):
            self.addChargeState(name=name, Z=Z, Z0=Z0, data=data[:,Z0,:])


    def __getitem__(self, Z0):
        """
        Returns the 'IonState' object corresponding to the
        specified charge state.
        """
        return self.ionstates[Z0]


    def __repr__(self):
        """
        Convert this IonSpecies object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this IonSpecies object to a string.
        """
        s = "Ion species {} (Z={}) of size NT x NR = {} x {}\n".format(self.name, self.Z, self.grid.t.size, self.grid.r.size)
        s += 't=final {}'.format(self.getDensity(t=-1))

        return s


    def addChargeState(self, name, Z, Z0, data):
        """
        Adds a new IonState object to the list of ion charge states.
        """
        self.ionstates.append(IonState(name=name, Z=Z, Z0=Z0, data=data, grid=self.grid, output=self.output))


    def getCharge(self): return self.Z


    def getName(self): return self.name


    def getDensity(self, t=-1):
        """
        Returns the total radial density (summed over all charge
        states) in the given time step.
        """
        n = np.zeros((self.grid.r.size,))

        for ion in self.ionstates:
            n += ion.get(t=t)

        return n


