# Base class for kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np
import scipy

from . OutputException import OutputException
from . UnknownQuantity import UnknownQuantity

from .. Settings.MomentumGrid import TYPE_PXI, TYPE_PPARPPERP


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
        return '({}) Kinetic quantity of size NT x NR x NP2 x NP1 = {} x {} x {} x {}\n:: {}\n:: Evolved using: {}\n{}'.format(self.name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3], self.description, self.description_eqn, self.data)


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def angleAveraged(self, t=None, r=None, moment='distribution'):
        """
        Returns the angle-averaged distribution function. Depending on
        the input parameters, the whole or only some parts of the spatiotemporal
        distribution can be angle-averaged.

        This method can only be applied to distributions defined on p/xi
        momentum grids.
        """
        if self.momentumgrid is None or self.momentumgrid.type != TYPE_PXI:
            raise OutputException("The angle average can only be calculated on p/xi grids.")

        if t is None: t = slice(None)
        if r is None: r = slice(None)
        
        data = self.data[t,r,:]

        if type(moment) == str:
            if moment == 'distribution': pass
            elif moment == 'density':
                data = data * self.momentumgrid.Vprime_VpVol[r,:]
            elif moment == 'current':
                if self.momentumgrid.type == TYPE_PXI:
                    absXi = abs(self.momentumgrid.p2[None,:,None])
                else: 
                    absXi = abs(np.atan(self.momentumgrid.p2[None,:,None]/self.momentumgrid.p1[None,None,:]))
                isPassing = absXi > self.momentumgrid.rgrid.xi0TrappedBoundary[r,None,None]
                data = data * isPassing * self.momentumgrid.getVpar() * self.momentumgrid.Vprime_VpVol[r,:] * scipy.constants.e
        elif type(moment) == float or type(moment) == np.ndarray:
            data = data * moment * self.momentumgrid.Vprime_VpVol[r,:]
        else:
            raise OutputException("Invalid type of parameter 'moment'.")
            
        favg = np.sum(data * self.momentumgrid.DP2[r,:], axis=data.ndim-2) / np.pi

        return favg


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


    def plot(self, t=-1, r=0, ax=None, show=None, logarithmic=False, coordinates='spherical', **kwargs):
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

        if coordinates is None:
            cp = ax.contourf(self.momentumgrid.p1, self.momentumgrid.p2, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(self.momentumgrid.getP1TeXName())
            ax.set_ylabel(self.momentumgrid.getP2TeXName())
        # Accept 'spherical' or 'spherica' or 'spheric' or ... 's':
        elif coordinates == 'spherical'[:len(coordinates)]:
            cp = ax.contourf(self.momentumgrid.P, self.momentumgrid.XI, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'$p$')
            ax.set_ylabel(r'$\xi$')
        elif coordinates == 'cylindrical'[:len(coordinates)]:
            cp = ax.contourf(self.momentumgrid.PPAR, self.momentumgrid.PPERP, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'$p_\parallel$')
            ax.set_ylabel(r'$p_\perp$')
        else:
            raise OutputException("Unrecognized coordinate type: '{}'.".format(coordinates))

        cb = None
        if genax:
            cb = plt.colorbar(mappable=cp, ax=ax)

        if show:
            plt.show(block=False)

        return ax


