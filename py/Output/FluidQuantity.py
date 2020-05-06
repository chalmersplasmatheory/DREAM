# Base class for fluid (radius + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np

from Output.UnknownQuantity import UnknownQuantity


class FluidQuantity(UnknownQuantity):
    

    def __init__(self, name, data, grid):
        """
        Constructor.
        """
        super(FluidQuantity, self).__init__(name=name, data=data, grid=grid)


    def plot(self, ax=None):
        """
        Generate a contour plot of the spatiotemporal evolution
        of this quantity.

        RETURNS a matplotlib axis object.
        """
        pass


    def plotRadialProfile(self, t=-1, ax=None):
        """
        Plot the radial profile of this quantity at the specified
        time slice.

        t: Time index to plot.

        RETURNS a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

        ax.plot(self.grid.r, self.data[t,:])
        ax.set_xlabel(r'Radius $r/a$')
        ax.set_ylabel('{}'.format(self.getTeXName()))

        return ax


    def plotTimeProfile(self, r=0, ax=None):
        """
        Plot the temporal profile of this quantity at the specified
        radius.

        r: Radial index to plot evolution for.

        RETURNS a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

        ax.plot(self.grid.t, self.data[:,r])
        ax.set_xlabel(r'Time $t$')
        ax.set_ylabel('{}'.format(self.getTeXName()))

        return ax
        

