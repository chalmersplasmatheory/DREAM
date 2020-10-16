# ppar/pperp output momentum grid object


import numpy as np
import scipy.constants
from .MomentumGrid import MomentumGrid


class PparPperpGrid(MomentumGrid):
    

    def __init__(self, name, rgrid, data):
        """
        Constructor.

        name:  Grid name.
        rgrid: Parent 'Grid' object (representing radial grid).
        data:  Momentum grid data.
        """
        super(PXiGrid, self).__init__(name=name, rgrid=rgrid, data=data)

        self.ppar   = data['p1']
        self.pperp  = data['p2']
        self.dppar  = data['dp1']
        self.dpperp = data['dp2']

        self.PPAR, self.PPERP = np.meshgrid(self.ppar, self.pperp)
        self.P = np.sqrt(self.PPAR**2 + self.PPERP**2)
        self.XI = self.PPAR / self.P
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


