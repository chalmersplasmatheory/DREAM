# p/xi output momentum grid object


import numpy as np
import scipy.constants
from .MomentumGrid import MomentumGrid


class PXiGrid(MomentumGrid):
    

    def __init__(self, name, rgrid, data):
        """
        Constructor.

        name:  Grid name.
        rgrid: Parent 'Grid' object (representing radial grid).
        data:  Momentum grid data.
        """
        super(PXiGrid, self).__init__(name=name, rgrid=rgrid, data=data)

        self.p1name = 'p'
        self.p2name = 'xi'

        self.p  = data['p1']
        self.xi = data['p2']
        self.dp = data['dp1']
        self.dxi = data['dp2']

        self.P, self.XI = np.meshgrid(self.p, self.xi)
        self.PPAR = self.P*self.XI
        self.PPERP = self.P*np.sqrt(1-self.XI**2)
        self.GAMMA = np.sqrt(self.P**2 + 1)


    def getGamma(self):
        """
        Returns a meshgrid representing the relativistic factor on this
        2D momentum grid.
        """
        return self.GAMMA


    def getVpar(self):
        """
        Returns a meshgrid representing the parallel velocity on this
        2D momentum grid.
        """
        return scipy.constants.c * (self.PPAR/self.GAMMA)


