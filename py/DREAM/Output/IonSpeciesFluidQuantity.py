# Representation of a fluid quantity which has an ion index.

import matplotlib.pyplot as plt
from . FluidQuantity import FluidQuantity
from . UnknownQuantity import UnknownQuantity


class IonSpeciesFluidQuantity(UnknownQuantity):
    

    def __init__(self, name, data, attr, grid, output):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)

        self.attr = attr
        self.ions = output.ionmeta


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        s = self.__str__()

        if hasattr(self, 'description') and hasattr(self, 'description_eqn'):
            s += "\n:: {}\n:: Evolved using: {}\n".format(self.description, self.description_eqn)

        return s


    def __str__(self):
        """
        Convert this object to a string.
        """
        if self.data.ndim == 3:
            nt, ni, nr = self.data.shape
        else:
            ni, nt, nr = 1, *self.data.shape
        s = '({}) Ion species fluid quantity of size NI x NT x NR = {} x {} x {}\n'.format(self.name, ni, nt, nr)
        for i in range(len(self.ions.Z)):
            s += "  {:2s} (Z = {:3d})\n".format(*self.ions[i])

        return s


    def __getitem__(self, name):
        """
        Direct access to data.
        """
        idx = self.ions.getIndex(name)

        if self.data.ndim == 3:
            data = self.data[:,idx,:]
        else:
            data = self.data[:]

        return FluidQuantity(name='{}_{}'.format(self.name, name), data=data, grid=self.grid, output=self.output, attr=self.attr)

    
    def dumps(self, ion=None, r=None, t=None):
        """
        Print the data in this quantity.
        """
        return self.get(ion=ion, r=r, t=t).__str__()


    def get(self, ion=None, r=None, t=None):
        """
        Returns data for the specified ion, or in the specified time
        interval or radial point. If none of the indices are given, returns
        the full evolution of the quantity.
        """
        sion = ion if ion is not None else slice(None)
        sr = r if r is not None else slice(None)
        st = t if t is not None else slice(None)

        if self.data.ndim == 3:
            return self.data[st,sion,sr]
        else:
            return self.data[st,sr]


    def plot(self, ion=None, ax=None, show=None, r=None, t=None, *args, **kwargs):
        """
        Plot data for all members of this IonSpeciesFluidQuantity in the
        same figure.
        """
        # Prevent trying to plot multiple 2D plots in the same window...
        if ion is None:
            if (r is None and self.grid.r.size != 1) and (t is None):
                raise Exception('Cannot plot ion temperature for all ions simultaneously when nr > 1.')

        if ion is not None:
            q = self[ion]
            ax = q.plot(ax=ax, show=show, r=r, t=t, *args, **kwargs)
        else:
            for i in self.ions.names:
                q = self[i]

                ax = q.plot(ax=ax, show=show, r=r, t=t, label=i, *args, **kwargs)

            plt.legend()

        return ax


