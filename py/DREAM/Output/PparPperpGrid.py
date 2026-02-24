# ppar/pperp output momentum grid object


import numpy as np
import scipy.constants
from .MomentumGrid import MomentumGrid


class PparPperpGrid(MomentumGrid):
    P1_NAME = 'ppar'
    P2_NAME = 'pperp'
    P1_TEX_NAME = r'$p_\parallel$'
    P2_TEX_NAME = r'$p_\perp$'

    def __init__(self, name, rgrid, data):
        """
        Constructor.

        name:  Grid name.
        rgrid: Parent 'Grid' object (representing radial grid).
        data:  Momentum grid data.
        """
        super(PparPperpGrid, self).__init__(name=name, rgrid=rgrid, data=data)

        self.ppar = self.p1
        self.pperp = self.p2
        self.ppar_f = self.p1_f
        self.pperp_f = self.p2_f
        self.dppar = self.dp1
        self.dpperp = self.dp2

        self._PPAR, self._PPERP = np.meshgrid(self.ppar[:], self.pperp[:])
        self._P = np.sqrt(self._PPAR**2 + self._PPERP**2)
        self._XI = self._PPAR / self._P
        self._GAMMA = np.sqrt(self._P**2 + 1)

        self._PPAR_f, self._PPERP_f = np.meshgrid(self.ppar_f[:], self.pperp_f[:])
        self._P_f = np.sqrt(self._PPAR_f**2 + self._PPERP_f**2)
        self._XI_f = self._PPAR_f / self._P_f

        self._VprimeCylindrical = np.copy(self.Vprime_VpVol)
        self._VprimeSpherical = self.Vprime_VpVol * self._P**2 / self._PPERP

    def getVpar(self):
        """
        Returns a meshgrid representing the parallel velocity on this
        2D momentum grid.
        """
        return scipy.constants.c * (self.PPAR/self.GAMMA)


