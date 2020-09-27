# Object for handling the ion density 'n_i'

import numpy as np
from .IonSpecies import IonSpecies
from .OutputException import OutputException
from .UnknownQuantity import UnknownQuantity


class IonHandler(UnknownQuantity):

    
    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super(IonHandler, self).__init__(name=name, data=None, attr=attr, grid=grid, output=output)

        self.meta = output.ionmeta
        self.ions = list()

        # Verify that data is consistent
        nZ0 = sum([Z+1 for Z in self.meta.Z])
        if nZ0 != data.shape[1]:
            raise OutputException("Inconsistent size of data array. Expected nZ0 = {}, but it was {}."
                .format(nZ0, data.shape[1]))

        # Add ions
        iidx = 0
        for name, Z in zip(self.meta.names, self.meta.Z):
            self.addIon(name=name, Z=Z, data=data[:,iidx:(iidx+Z+1),:], attr=attr)
            iidx += Z+1


    def __getitem__(self, i):
        """
        If i is a string, retrieves an 'IonSpecies' object by name.
        Otherwise, if i is an integer, returns the ion by index.

        If also Z0 is specified (and is an integer), returns an
        'IonState' object corresponding to the selected charge state.
        """
        idx, Z0 = None, None

        ion = None
        if type(i) == str:
            ion = self.getIonByName(i)
        else:
            if len(i) == 2:
                idx = i[0]
                Z0 = i[1]
            else:
                idx = i

            ion = self.ions[idx]

        if Z0 is not None:
            return ion[Z0]
        else:
            return ion


    def __repr__(self):
        """
        Convert this IonHandler to an "official" string representation.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this IonHandler object to a string.
        """
        # TODO also show relative abundance
        s = "{} ion species\n".format(len(self.ions))
        for ion in self.ions:
            s += "   {:2s} (Z = {:3d})  {:.3e} particles\n".format(ion.getName(), ion.getCharge(), ion.getParticleNumber())
        
        return s


    def addIon(self, name, Z, data, attr=list()):
        """
        Adds a new ion to the list of ions.
        """
        self.ions.append(IonSpecies(name=name, Z=Z, data=data, grid=self.grid, output=self.output, attr=attr))
    
    
    def getIonByName(self, name):
        """
        Returns the ion with the specified name.
        """
        for ion in self.ions:
            if ion.getName() == name:
                return ion

        raise KeyError("No ion named '{}' found in the output.".format(name))

