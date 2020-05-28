# Base class for "other" kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np

from . OutputException import OutputException
from . OtherQuantity import OtherQuantity


class OtherKineticQuantity(OtherQuantity):
    

    def __init__(self, name, data, grid, output, momentumgrid=None):
        """
        Constructor.
        """
        super(OtherKineticQuantity, self).__init__(name=name, data=data, grid=grid, output=output)

        self.momentumgrid = momentumgrid


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        #s = self.__str__() 
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Other kinetic quantity of size NT x NR x NP2 x NP1 = {} x {} x {} x {}'.format(self.name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3])


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


    def getCoordinateGrids(self):
        """
        Returns the appropriate coordinate grid vectors
        for this other quantity. Which grid vectors to
        return is determined based on the name of the quantity.
        If the name ends with

          '_fr': Returns radial flux grid vector + distribution momentum vectors
          '_f1': Returns p1 flux grid vector + distribution r and p2 vectors
          '_f2': Returns p2 flux grid vector + distribution r and p1 vectors
        """
        suff = self.name[-3:]
        r  = self.grid.r
        p1 = self.momentumgrid.p1
        p2 = self.momentumgrid.p2

        if   suff == '_fr': r  = self.grid.r_f
        elif suff == '_f1': p1 = self.momentumgrid.p1_f
        elif suff == '_f2': p2 = self.momentumgrid.p2_f

        return r, p1, p2


    def plot2D(self, t=-1, r=0, ax=None, show=None):
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

        data = self.data[t,r,:]

        if data.ndim != 2:
            raise OutputException("Data dimensionality is too high. Unable to visualize kinetic quantity.")

        _, p1, p2 = self.getCoordinateGrids()

        cp = ax.contourf(p1, p2, data, cmap='GeriMap')
        ax.set_xlabel(self.momentumgrid.getP1TeXName())
        ax.set_ylabel(self.momentumgrid.getP2TeXName())

        cb = None
        if genax:
            cb = plt.colorbar(mappable=cp, ax=ax)

        if show:
            plt.show(block=False)

        return ax


