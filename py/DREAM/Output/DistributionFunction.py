# Distribution function data type
#

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants

from . KineticQuantity import KineticQuantity
from . OutputException import OutputException
from .. import GeriMap
from .. Settings.MomentumGrid import TYPE_PXI, TYPE_PPARPPERP


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
            if self.momentumgrid.type == TYPE_PXI:
                p1name, p2name = 'P', 'XI'
            elif self.momentumgrid.type == TYPE_PPARPPERP:
                p1name, p2name = 'PAR', 'PERP'

        return '({}) Kinetic quantity of size NT x NR x N{} x N{} = {} x {} x {} x {}'.format(self.name, p2name, p1name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3])


    #########################################
    # INTEGRALS OF THE DISTRIBUTION FUNCTION
    #########################################
    def angleAveraged(self, t=None, r=None):
        """
        Returns the angle-averaged distribution function. Depending on
        the input parameters, the whole or only some parts of the spatiotemporal
        distribution can be angle-averaged.

        This method can only be applied to distributions defined on p/xi
        momentum grids.
        """
        if self.momentumgrid is None or self.momentumgrid.type != TYPE_PXI:
            raise OutputException("The distribution angle average can only be calculated on p/xi grids.")

        data = self.data[t,r,:]
        favg = np.sum(data, axis=data.ndim-2) / np.pi

        return favg


    def currentDensity(self, t=None, r=None):
        """
        Calculates the current density carried by the electrons of
        this distribution function.
        """
        if t is None:
            t = range(len(self.grid.t))
        if r is None:
            r = range(len(self.grid.r))

        if np.isscalar(t):
            t = np.asarray([t])
        if np.isscalar(r):
            r = np.asarray([r])

        Vpar = self.momentumgrid.getVpar()

        j = []
        for iT in range(len(t)):
            jr = []
            for iR in range(len(r)):
                jr.append(self.momentumgrid.integrate2D(self.data[t[iT],r[iR],:] * Vpar)[0])

            j.append(jr)

        j = np.asarray(j) * scipy.constants.e

        return j


    def density(self, t=None, r=None):
        """
        Calculates the total density of this distribution function.
        """
        if t is None:
            t = range(len(self.grid.t))
        if r is None:
            r = range(len(self.grid.r))

        if np.isscalar(t):
            t = np.asarray([t])
        if np.isscalar(r):
            r = np.asarray([r])

        n = []
        for iT in range(len(t)):
            nr = []
            for iR in range(len(r)):
                nr.append(self.momentumgrid.integrate2D(self.data[t[iT],r[iR],:])[0])

            n.append(nr)

        n = np.asarray(n)

        return n


    def plasmaCurrent(self, t=None):
        """
        Calculates the total plasma current carried by the electrons of
        this distribution function.
        """
        j = self.currentDensity(t=t)
        return self.grid.integrate(j)
        #if t is None:
        #    return self.momentumgrid.integrate3D(self.data, 
        #else:
        #    return self.momentumgrid.integrate3D(self.data[t,:]) 


    ##########################################
    # PLOTTING ROUTINES
    ##########################################
    def plot(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Alias for 'semilogy()' henceforth.
        """
        return self.semilogy(t=t, r=r, p2=p2, ax=ax, show=show)


    def plot2D(self, t=-1, r=0, ax=None, show=None, logarithmic=True):
        """
        Make a contour plot of this quantity.
        """
        return super(DistributionFunction, self).plot(t=t, r=r, ax=ax, show=show, logarithmic=logarithmic)


    def semilog(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Alias for 'semilogy()'.
        """
        return self.semilogy(t=t, r=r, p2=p2, ax=ax, show=show)


    def semilogy(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Plot this distribution function on a semilogarithmic scale.
        If 'p2' is None, the distribution function is first angle-averaged.
        Otherwise, 'p2' is interpreted as an index into the distribution
        function (second momentum dimension, i.e. xi or pperp).
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
            favg = self.data[slice(t),slice(r),slice(p2),:]

        #if favg.ndim != 1:
        #    raise OutputException("Data dimensionality is too high. Unable to visualize distribution function.")
        if favg.ndim > 1: ndim = np.prod(favg.shape[:-1])
        else: ndim = 1

        favg = np.reshape(favg, (ndim, favg.shape[-1]))

        colors = GeriMap.get(N=ndim+1)
        lbls = []
        for i in range(0, ndim):
            ax.semilogy(self.momentumgrid.p1, favg[i,:], color=colors(i/(ndim+1)))

            if np.isscalar(t) and np.isscalar(r): continue
            elif np.isscalar(r):
                tval, unit = self.grid.getTimeAndUnit(t[i])
                lbls.append(r'$t = {:.3f}\,\mathrm{{{}}}$'.format(tval, unit))
            elif np.isscalar(t):
                lbls.append(r'$r = {:.3f}\,\mathrm{{m}}$'.format(self.grid.r[r[i]]))
            else:
                it, ir = np.unravel_index(i, (len(t), len(r)))
                tval, unit = self.grid.getTimeAndUnit(it)
                lbls.append(r'$t = {:.3f}\,\mathrm{{{}}}, r = {:.3f}\,\mathrm{{m}}$'.format(tval, unit, self.grid.r[ir]))

        ax.set_xlabel(self.momentumgrid.getP1TeXName())
        ax.set_ylabel(self.getTeXName())

        fmax = np.amax(favg)
        ax.set_xlim([self.momentumgrid.p1[0], self.momentumgrid.p1[-1]])
        ax.set_ylim(np.array([1e-30, 10]) * fmax)

        if len(lbls) > 0:
            ax.legend(lbls)

        if show:
            plt.show(block=False)
        
        return ax


