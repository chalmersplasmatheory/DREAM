# Object for handling the ion density 'n_i'

import matplotlib.pyplot as plt
import numpy as np

from .IonSpecies import IonSpecies
from .OutputException import OutputException
from .UnknownQuantity import UnknownQuantity


class IonHandler(UnknownQuantity):

    
    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super(IonHandler, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output)

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


    def getIonOffset(self, name, Z0=0):
        """
        Returns the array index of the named ion species (and charge state)
        which can be used to index a vector with all ion species located one
        after another.
        """
        return self.meta.getIonOffset(name=name, Z0=Z0)


    def ionNameToIndex(self, name):
        """
        Returns the index of the named ion species.
        """
        for i in range(len(self.ions)):
            if self.ions[i].getName() == name:
                return i

        raise KeyError("No ion named '{}' found in the output.".format(name))


    def plot(self, t=-1, ax=None, show=None, **kwargs):
        """
        Visualize the ion charge state densities.
        """
        self.histogram(t=t, ax=ax, show=show, **kwargs)


    def histogram(self, t=-1, ax=None, show=None, **kwargs):
        """
        Create a histogram of the ion charge state densities.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        x = [0]
        y = [0]
        labels = ['']
        for i in range(len(self.ions)):
            ion = self.ions[i]
            for state in ion.ionstates:
                x.append(x[-1]+1)
                y.append(state.integral(t=t))
                labels.append(state.getRomanName())
        
        x = x[1:]
        y = y[1:]
        labels = labels[1:]

        #ax.hist(y, weights=x)
        ax.bar(x, y, tick_label=labels, **kwargs)
        lbls = ax.get_xticklabels()
        plt.setp(lbls, rotation=45)

        if show:
            plt.show(block=False)

        return ax


