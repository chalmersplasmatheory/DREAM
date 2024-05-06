# Distribution function data type
#

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants

from . import Bekefi
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

    def currentDensity(self, t=None, r=None):
        """
        Calculates the current density carried by the electrons of
        this distribution function.
        """
#        Vpar = self.momentumgrid.getBounceAveragedVpar()
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


    def plasmaCurrent(self, t=None):
        """
        Calculates the total plasma current carried by the electrons of
        this distribution function.
        """
        j = self.currentDensity(t=t)
        return self.grid.integrate(j)


    def pressure(self, t=None, r=None):
        """
        Evaluates the parallel and perpendicular pressure moments of the
        distribution function according to

          P_||    = mc^2 <p_||^2 / sqrt(1+p^2)>
          P_\perp = mc^2 <p_\perp^2 / (2*sqrt(1+p^2))>

        where <...> denotes an integral over all of momentum space. The
        pressure is returned in units of Pa.

        :returns: Tuple consisting of parallel and perpendicular pressures.
        """
        p, xi  = self.momentumgrid.p1[:], self.momentumgrid.p2[:]
        P, XI  = np.meshgrid(p, xi)
        PPAR2  = (P*XI)**2
        PPERP2 = P**2*(1-XI**2)

        GAMMA = self.momentumgrid.getGamma()
        mc2 = scipy.constants.m_e * scipy.constants.c**2

        Ppar = self.moment(PPAR2/GAMMA, t=t, r=r) * mc2
        Pperp = self.moment(PPERP2/GAMMA, t=t, r=r) * 0.5*mc2

        return Ppar, Pperp


    def partialDensity(self, t=None, r=None, pc=None):
        """
        Calculate the density of electrons contained above a momentum ``pc``. If
        a value for ``pc`` is not explicitly given, it is taken as either the
        critical momentum ``other.fluid.pCrit`` if available, or as the
        approximation ``pc = 1/sqrt(E/Ec-1)``.

        The parameter ``pc`` may be a scalar or a fluid quantity, varying in
        time and space.
        """
        if t is None:
            t = list(range(len(self.time)))
        elif np.isscalar(t):
            t = np.array([t])

        if r is None:
            r = list(range(len(self.grid.r)))
        elif np.isscalar(t):
            r = np.array([r])

        if pc is None:
            if 'fluid' in self.output.other:
                if 'pCrit' in self.output.other.fluid:
                    pc = self.output.other.fluid.pCrit[:]
                elif 'Ecfree' in self.output.other.fluid:
                    Ec = self.output.other.fluid.Ecfree[:]
                    E = self.output.eqsys.E_field[1:,:]
                    pc = 1/np.sqrt(E/Ec - 1)

            if pc is None:
                T = self.output.eqsys.T_cold[:]
                n = self.output.eqsys.n_cold[:]
                E = self.output.eqsys.E_field[:]

                Ec = Formulas.Ec(T=T, n=n)
                pc = 1/np.sqrt(E/Ec-1)

            if pc.shape[0] == self.grid.t.size-1:
                ppc = np.zeros((self.grid.t.size, pc.shape[1]))
                ppc[1:,:] = pc[:]
                ppc[0,:]  = np.inf
                pc = ppc
        elif np.isscalar(pc):
            pc = pc*np.ones((self.time.size, self.grid.r.size))

        p = self.momentumgrid.p1
        data = np.zeros((len(t), len(r)))
        for i in range(len(t)):
            for j in range(len(r)):
                F = self.data
                mask = np.zeros((F.shape[2], F.shape[3]))
                mask[:,np.where(p >= pc[t[i],r[j]])] = 1
                data[i,j] = self.moment(mask, t=t[i], r=r[j])

        return data


    def runawayRate(self, t=None, r=None, pc=None):
        """
        Calculate the runaway rate, with the given definition of critical
        momentum, as a moment of the distribution function.
        """
        n = self.partialDensity(pc=pc)
        dt = np.diff(self.grid.t)
        dndt = (np.diff(n, axis=0).T / dt).T

        return dndt


    def synchrotron(self, model='spectrum', B=3.1, wavelength=700e-9, t=None, r=None):
        """
        Returns the synchrotron radiation emitted by this distribution function
        as a :py:class:`DREAM.Output.KineticQuantity.KineticQuantity`.

        :param str model:        Model to use for synchrotron moment. Either 'spectrum' (for synchrotron spectrum at specified wavelength and magnetic field) or 'total' (for total emitted power).
        :param float B:          Magnetic field strength to use with 'spectrum' model.
        :param float wavelength: Wavelength to use with 'spectrum' model (in meters).
        """
        if t is None:
            t = range(len(self.time))
        elif np.isscalar(t):
            t = np.array([t])

        if r is None:
            r = range(len(self.grid.r))
        elif np.isscalar(r):
            r = np.array([r])

        data = None
        if model == 'total':
            data = np.zeros((len(t), len(r), self.data.shape[-2], self.data.shape[-1]))
            pperp2 = self.momentumgrid.PPERP**2
            m2c2   = (scipy.constants.m_e * scipy.constants.c)**2
            for i in range(len(t)):
                for j in range(len(r)):
                    data[i,j,:,:] = pperp2 * self.data[t[i],r[j],:,:] * self.momentumgrid.Vprime_VpVol[r[j],:,:]
        elif model == 'spectrum':
            S = []
            W = Bekefi.synchrotron(self.momentumgrid.P, self.momentumgrid.XI, wavelength, B)
            data = W * self.data[t,r,:] * self.momentumgrid.Vprime_VpVol[r,:]
        else:
            raise OutputException("Unrecognized model for calculating synchrotron moment with: '{}'.".format(model))
        return KineticQuantity('synchrotron({})'.format(self.name), data=data, grid=self.grid, output=self.output, momentumgrid=self.momentumgrid, attr={'description': 'Synchrotron moment of {}'.format(self.name), 'equation': 'synchrotron({})'.format(self.name)})


    ##########################################
    # PLOTTING ROUTINES
    ##########################################
    def plot(self, t=-1, r=0, moment='distribution', p2=None, ax=None, show=None, logy=True, **kwargs):
        """
        Alias for :py:meth:`semilogy` henceforth.
        """
        v = self.semilogy(t=t, r=r, moment=moment, p2=p2, ax=ax, show=show, **kwargs)

        if logy:
            v.set_yscale('log')
        else:
            v.set_yscale('linear')

        return v


    def plot2D(self, t=-1, r=0, ax=None, show=None, logarithmic=True, coordinates=None, **kwargs):
        """
        Make a contour plot of this quantity.

        :param int t:            Time index to plot.
        :param int r:            Radial index to plot.
        :param ax:               Matplotlib axes object to use for plotting.
        :param bool show:        If ``True``, or ``None`` and ``ax`` is NOT provided, calls ``matplotlib.pyplot.show()`` just before returning.
        :param bool logarithmic: If ``True``, uses a logarithmic colour scale.
        :param str coordinates:  Name of coordinates to use (either 'spherical' (p/xi) or 'cylindrical' (ppar/pperp)).
        :param kwargs:           Keyword arguments passed on to matplotlib.contourf().
        """
        return super(DistributionFunction, self).plot(t=t, r=r, ax=ax, show=show, logarithmic=logarithmic, coordinates=coordinates, **kwargs)


    def semilog(self, t=-1, r=0, p2=None, ax=None, show=None, **kwargs):
        """
        Alias for :py:meth:`semilogy`.
        """
        return self.semilogy(t=t, r=r, p2=p2, ax=ax, show=show, **kwargs)


    def semilogy(self, t=-1, r=0, moment='distribution', p2=None, ax=None, show=None, **kwargs):
        """
        Plot this distribution function on a semilogarithmic scale.
        If 'p2' is None, the distribution function is first angle-averaged.
        Otherwise, 'p2' is interpreted as an index into the distribution
        function (second momentum dimension, i.e. xi or pperp).

        :param int t: Integer, or list of integers, specifying the time indices of the distributions to plot.
        :param int r: Integer, or list of integers, specifying the radial indices of the distributions to plot.
        :param str moment: String (or array) speciyfing the angle averaged moment of the distribution to plot. See :py:meth:`DREAM.Output.KineticQuantity.KineticQuantity.angleAveraged` for possible values.
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
        p = self.momentumgrid.p1
        for i in range(0, ndim):
            if 'color' not in kwargs:
                ax.semilogy(p, favg[i,:], color=colors(i/(ndim+1)), **kwargs)
            else:
                ax.semilogy(p, favg[i,:], **kwargs)

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


