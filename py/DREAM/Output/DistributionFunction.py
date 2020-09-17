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
    

    def __init__(self, name, data, grid, output, momentumgrid=None, attr=list()):
        """
        Constructor.
        """
        super(DistributionFunction, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output, momentumgrid=momentumgrid)


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

        return '({}) Kinetic quantity of size NT x NR x N{} x N{} = {} x {} x {} x {}\n:: {}\n:: Evolved using: {}'.format(self.name, p2name, p1name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3], self.description, self.description_eqn)


    #########################################
    # INTEGRALS OF THE DISTRIBUTION FUNCTION
    #########################################
    def angleAveraged(self, t=None, r=None, moment='distribution'):
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

        if type(moment) == str:
            if moment == 'distribution': pass
            elif moment == 'density':
                data = data * self.momentumgrid.Vprime_VpVol
            elif moment == 'current':
                data = data * self.momentumgrid.getVpar() * self.momentumgrid.Vprime_VpVol * scipy.constants.e
        elif type(moment) == float or type(moment) == np.ndarray:
            data = data * moment * self.momentumgrid.Vprime_VpVol
        else:
            raise OutputException("Invalid type of parameter 'moment'.")
            
        favg = np.sum(data * self.momentumgrid.DP2[r,:], axis=data.ndim-2) / np.pi

        return favg


    def currentDensity(self, t=None, r=None):
        """
        Calculates the current density carried by the electrons of
        this distribution function.
        """
        Vpar = self.momentumgrid.getVpar()
        return self.moment(Vpar, t=t, r=r) * scipy.constants.e


    def density(self, t=None, r=None):
        """
        Calculates the total density of this distribution function.
        """
        return self.moment(1, t=t, r=r)


    def kineticEnergy(self, t=None, r=None):
        """
        Calculates the kinetic energy contained in the distribution function.
        (energy in Joule).
        """
        gamma1 = self.momentumgrid.getGamma()-1
        mc2    = scipy.constants.m_e * scipy.constants.c**2
        return self.moment(gamma1, t=t, r=r) * mc2


    def moment(self, weight, t=None, r=None):
        """
        Evaluate a moment of this distribution function with the
        given weighting factor.
        """
        if t is None:
            t = range(len(self.grid.t))
        if r is None:
            r = range(len(self.grid.r))

        if np.isscalar(t):
            t = np.asarray([t])
        if np.isscalar(r):
            r = np.asarray([r])

        q = []
        for iT in range(len(t)):
            qr = []
            for iR in range(len(r)):
                qr.append(self.momentumgrid.integrate2D(self.data[t[iT],r[iR],:] * weight)[0])

            q.append(qr)

        q = np.asarray(q)

        return q


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
    def plot(self, t=-1, r=0, moment='distribution', p2=None, ax=None, show=None, logy=True):
        """
        Alias for 'semilogy()' henceforth.
        """
        v = self.semilogy(t=t, r=r, moment=moment, p2=p2, ax=ax, show=show)

        if logy:
            v.set_yscale('log')
        else:
            v.set_yscale('linear')

        return v


    def plot2D(self, t=-1, r=0, ax=None, show=None, logarithmic=True, **kwargs):
        """
        Make a contour plot of this quantity.
        """
        return super(DistributionFunction, self).plot(t=t, r=r, ax=ax, show=show, logarithmic=logarithmic, **kwargs)


    def semilog(self, t=-1, r=0, p2=None, ax=None, show=None):
        """
        Alias for 'semilogy()'.
        """
        return self.semilogy(t=t, r=r, p2=p2, ax=ax, show=show)


    def semilogy(self, t=-1, r=0, moment='distribution', p2=None, ax=None, show=None):
        """
        Plot this distribution function on a semilogarithmic scale.
        If 'p2' is None, the distribution function is first angle-averaged.
        Otherwise, 'p2' is interpreted as an index into the distribution
        function (second momentum dimension, i.e. xi or pperp).

        :param int t: Integer, or list of integers, specifying the time indices of the distributions to plot.
        :param int r: Integer, or list of integers, specifying the radial indices of the distributions to plot.
        :param str moment: String (or array) speciyfing the angle averaged moment of the distribution to plot. See :py:method:`angleAveraged` for possible values.
        :param int p2: Index into second momentum parameter (xi or pperp) of distribution to plot.
        :param matplotlib.Axes ax: Axes object to draw on.
        :param bool show: If ``True``, show the figure after plotting.
        """
        genax = ax is None

        # Generate matplotlib axes
        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        # Retrieve data to plot
        if p2 is None:
            favg = self.angleAveraged(t=t, r=r, moment=moment)
        else:
            favg = self.data[t,r,p2,:]

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

        if type(moment) == str:
            n = self.getTeXName()
            if moment == 'distribution':
                ax.set_ylabel(r'$\langle${}$\rangle$'.format(n))
            elif moment == 'density':
                ax.set_ylabel(r"$\langle(\mathcal{{V}}'/V')${}$\rangle$".format(n))
            elif moment == 'current':
                ax.set_ylabel(r"$\langle ev_\parallel(\mathcal{{V}}'/V')${}$\rangle$".format(n))
            else:
                ax.set_ylabel(r"Moment of $(\mathcal{{V}}'/V')${}".format(n))
        else:
            ax.set_ylabel(r"Moment of $(\mathcal{{V}}'/V')${}".format(n))

        fmax = np.amax(favg)
        ax.set_xlim([self.momentumgrid.p1[0], self.momentumgrid.p1[-1]])
        ax.set_ylim(np.array([1e-30, 10]) * fmax)

        if len(lbls) > 0:
            ax.legend(lbls)

        if show:
            plt.show(block=False)
        
        return ax


