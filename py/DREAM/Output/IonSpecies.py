# Object representing a single ion species

import matplotlib.pyplot as plt
import numpy as np
from .IonState import IonState
from .OutputException import OutputException
from .FluidQuantity import FluidQuantity


class IonSpecies:
    

    def __init__(self, name, Z, data, grid, output, attr=list()):
        """
        Constructor.

        name: Name of ion species.
        Z:    Ion atomic charge.
        data: Ion density data (size nt x nZ0 x nr, where nZ0 = Z+1
              is the number of charge states for this ion).
        """
        self.name = name
        self.Z    = Z
        self.grid = grid
        self.output = output
        self.ionstates = list()
        self.data = data

        nt = len(self.grid.t)
        nr = len(self.grid.r)

        if data.shape != (nt, Z+1, nr):
            raise OutputException("Invalid dimensions of data. Expected {}x{}x{}, but found {}x{}x{}."
                .format(nt, Z+1, nr, data.shape[0], data.shape[1], data.shape[2]))

        for Z0 in range(0, Z+1):
            self.addChargeState(name=name, Z=Z, Z0=Z0, data=data[:,Z0,:], attr=attr)


    def __getitem__(self, Z0):
        """
        Returns the 'IonState' object corresponding to the
        specified charge state.
        """
        return self.ionstates[Z0]


    def __repr__(self):
        """
        Convert this IonSpecies object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this IonSpecies object to a string.
        """
        s = "Ion species {} (Z={}) of size NT x NR = {} x {}\n".format(self.name, self.Z, self.grid.t.size, self.grid.r.size)
        #s += 't=final {}'.format(self.getDensity(t=-1))
        for Z0 in range(self.Z+1):
            s += "  Z0={:2d}:  {:1.3e} particles\n".format(Z0, self.ionstates[Z0].integral(t=-1))

        return s


    def addChargeState(self, name, Z, Z0, data, attr=list()):
        """
        Adds a new IonState object to the list of ion charge states.
        """
        self.ionstates.append(IonState(name=name, Z=Z, Z0=Z0, data=data, grid=self.grid, output=self.output, attr=attr))


    def getCharge(self): return self.Z


    def getName(self): return self.name


    def getDensity(self, t=-1):
        """
        Returns the total radial density (summed over all charge
        states) in the given time step.
        """
        #n = np.zeros((self.grid.r.size,))
        n = None

        for ion in self.ionstates:
            if n is None:
                n = ion.get(t=t)
            else:
                n += ion.get(t=t)

        return n


    def getParticleNumber(self, t=-1):
        """
        Returns the number of particles of this ion species
        in the given time step.
        """
        return self.grid.integrate(self.getDensity(t=t))


    def plot(self, t=None, r=None, Z0=None, ax=None, show=None):
        """
        Plots the ion charge state densities.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        if t is None: t = slice(None)
        if Z0 is None: Z0 = slice(None)

        legs = []
        states = self.ionstates[Z0]
        for state in states:
            data = None
            # If 'r' is None, we integrate over all radii
            # to get the total number of particles...
            if r is None:
                data = state.integral(t=t)
            else:
                data = state[t,r]

            ax.plot(self.grid.t[t], data)

            legs.append('$Z_0 = {}$'.format(state.Z0))

        if len(legs) > 0:
            ax.legend(legs)

        if show:
            plt.show(block=False)

        return ax
        
    def plotSum(self, Z0 = None, integrate = False, **kwargs):
        """
        Plots the spatio-temporal evolution of the sum of the specified charge states of this ion species
        Z0:         list of charge states to be summed and plotted
        integrate:  if 'True', plot the volume integral of the specified charge states
        """
        if Z0 is None: Z0 = slice(None)
        states = self.ionstates[Z0]
        data = None
        for state in states:
            if data is None:
                data = state.data
            else:
                data = data + state.data
                
        _IonSum = FluidQuantity(name = self.name, data=data, attr=list(), grid=self.grid, output=self.output)
        if integrate:
            return _IonSum.plotIntegral(**kwargs)
        else:
            return _IonSum.plot(**kwargs)
        
        
        
        
