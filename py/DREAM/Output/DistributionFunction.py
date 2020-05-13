# Distribution function data type
#

import numpy as np
import matplotlib.pyplot as plt

from . KineticQuantity import KineticQuantity
from . OutputException import OutputException
from .. Settings.MomentumGrid import MOMENTUMGRID_TYPE_PXI, MOMENTUMGRID_TYPE_PPARPPERP


class DistributionFunction(KineticQuantity):
    

    def __init__(self, name, data, grid, output, momentumgrid=None):
        """
        Constructor.
        """
        super(DistributionFunction, self).__init__(name=name, data=data, grid=grid, output=output, momentumgrid=momentumgrid)


    def __str__(self):
        """
        Convert this object to a string.
        """
        p1name = 'P1'
        p2name = 'P2'

        if self.momentumgrid is not None:
            if self.momentumgrid.type == MOMENTUMGRID_TYPE_PXI:
                p1name, p2name = 'P', 'XI'
            elif self.momentumgrid.type == MOMENTUMGRID_TYPE_PPARPPERP:
                p1name, p2name = 'PAR', 'PERP'

        return '({}) Kinetic quantity of size NT x NR x N{} x N{} = {} x {} x {} x {}'.format(self.name, p2name, p1name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3])


    def angleAveraged(self, t=None, r=None):
        """
        Returns the angle-averaged distribution function. Depending on
        the input parameters, the whole or only some parts of the spatiotemporal
        distribution can be angle-averaged.

        This method can only be applied to distributions defined on p/xi
        momentum grids.
        """
        if self.momentumgrid is None or self.momentumgrid.type != MOMENTUMGRID_TYPE_PXI:
            raise OutputException("The distribution angle average can only be calculated on p/xi grids.")

        data = self.data[t,r,:,:]
        favg = np.sum(data, axis=data.ndim-2) / np.pi

        return favg


    def plot(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Alias for 'semilogy()' henceforth.
        """
        return self.semilogy(t=t, r=r, p2=p2, ax=ax, show=show)


    def plot2D(self, t=-1, r=0, ax=None, show=None):
        """
        Make a contour plot of this quantity.
        """
        return super(DistributionFunction, self).plot(t=t, r=r, ax=ax, show=show)


    def semilog(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Alias for 'semilogy()'.
        """
        return self.semilogy(t=t, r=r, p2=p2, ax=ax, show=show)


    def semilogy(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Plot this distribution function on a semilogarithmic scale.
        If 'p2' is None, the distribution function is first angle-averaged.
        """
        genax = ax is None

        # Generate matplotlib axes
        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        # Retrieve data to plot
        if p2 is None:
            favg = self.angleAveraged(t=t, r=r)
        else:
            favg = self.data[t,r,p2,:]

        if favg.ndim != 1:
            raise OutputException("Data dimensionality is too high. Unable to visualize distribution function.")

        cp = ax.semilogy(self.momentumgrid.p1, favg)
        ax.set_xlabel(self.momentumgrid.getP1TeXName())
        ax.set_ylabel(self.getTeXName())

        fmax = np.amax(favg)
        ax.set_xlim([self.momentumgrid.p1[0], self.momentumgrid.p1[-1]])
        ax.set_ylim(np.array([1e-30, 10]) * fmax)

        if show:
            plt.show(block=False)
        
        return ax


