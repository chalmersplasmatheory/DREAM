# Base class for kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np

from . OutputException import OutputException
from . UnknownQuantity import UnknownQuantity


class KineticQuantity(UnknownQuantity):
    

    def __init__(self, name, data, grid, output, momentumgrid=None, attr=list()):
        """
        Constructor.
        """
        super(KineticQuantity, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output)

        self.momentumgrid = momentumgrid


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Kinetic quantity of size NT x NR x NP2 x NP1 = {} x {} x {} x {}\n:: {}\n:: Evolved using: {}\n'.format(self.name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3], self.description, self.description_eqn)


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def get(self, t=None, r=None, p2=None, p1=None):
        """
        Returns data using the specified indexing. If an
        argument is 'None', this method will return all
        elements along that dimension.
        """
        sel = [slice(None)] * 4

        if t  is not None: sel[0] = t
        if r  is not None: sel[1] = r
        if p2 is not None: sel[2] = p2
        if p1 is not None: sel[3] = p1

        return self.data[tuple(sel)]


    def plot(self, t=-1, r=0, ax=None, show=None, logarithmic=False, **kwargs):
        """
        Plot this kinetic quantity.
        """
        if self.momentumgrid is None:
            raise OutputException("Unable to plot kinetic quantity as its momentum grid has not been specified.")

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        data = None
        if logarithmic:
            data = np.log10(self.data[t,r,:])
        else:
            data = self.data[t,r,:]

        if data.ndim != 2:
            raise OutputException("Data dimensionality is too high. Unable to visualize kinetic quantity.")

        cp = ax.contourf(self.momentumgrid.p1, self.momentumgrid.p2, data, cmap='GeriMap', **kwargs)
        ax.set_xlabel(self.momentumgrid.getP1TeXName())
        ax.set_ylabel(self.momentumgrid.getP2TeXName())

        cb = None
        if genax:
            cb = plt.colorbar(mappable=cp, ax=ax)

        if show:
            plt.show(block=False)

        return ax


