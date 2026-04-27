# p/xi output momentum grid object


import numpy as np
import scipy.constants
from .MomentumGrid import MomentumGrid


class PXiGrid(MomentumGrid):
    P1_NAME = 'p'
    P2_NAME = 'xi'
    P1_TEX_NAME = r'$p$'
    P2_TEX_NAME = r'$\xi$'

    def __init__(self, name, rgrid, data):
        """
        Constructor.

        name:  Grid name.
        rgrid: Parent 'Grid' object (representing radial grid).
        data:  Momentum grid data.
        """
        super(PXiGrid, self).__init__(name=name, rgrid=rgrid, data=data)
        self.p  = self.p1
        self.xi = self.p2
        self.p_f = self.p1_f
        self.xi_f = self.p2_f
        self.dp = self.dp1
        self.dxi = self.dp2

        def ppar(p, xi):
            return p * xi
        def pperp(p, xi):
            sinSq = 1-xi**2
            sinSq[sinSq<0] = 0
            return p * np.sqrt(sinSq)

        self._P, self._XI = np.meshgrid(self.p[:], self.xi[:])
        self._PPAR = ppar(self._P, self._XI)
        self._PPERP = pperp(self._P, self._XI)
        self._GAMMA = np.sqrt(self._P**2 + 1)

        # flux grid corners grid, for use with e.g. pcolormesh
        self._P_f, self._XI_f = np.meshgrid(self.p_f[:], self.xi_f[:])
        self._PPAR_f = ppar(self._P_f, self._XI_f)
        self._PPERP_f = pperp(self._P_f, self._XI_f)

        # phase-space weights in various representations
        self._VprimeCylindrical = self.Vprime_VpVol * self._PPERP / self._P**2
        self._VprimeSpherical = np.copy(self.Vprime_VpVol)

    def getBounceAveragedVpar(self):
        """
        Returns a meshgrid representing the integrand that should weigh
        a function when carrying out the v_par moment of a quantity.
        
        It should be identical to the ``integrand`` produced by the 
        ``CurrentDensityFromDistributionFunction`` class in the DREAM kernel.        
        """
        c  = scipy.constants.speed_of_light
        integrand = np.zeros(self.Vprime.shape)

        # Load in data from file to speed up calculation
        Vprime_VpVol = self.Vprime_VpVol[:]
        p = self.p[:]
        xi = self.xi[:]
        xi_f = self.xi_f[:]
        xi0TrappedBoundary = self.rgrid.xi0TrappedBoundary[:]
        v = c * p / np.sqrt(1+p**2)

        # Calculate bounce averaged v||
        for ir in range(0, self.rgrid.r.size):
            xi0Trapped = xi0TrappedBoundary[ir]
            for j in range(0, xi.size):
                xi1 = xi_f[j]
                xi2 = xi_f[j+1]
                if(xi1>xi2):
                    xi_t = xi1
                    xi1 = xi2
                    xi2 = xi_t
                
                xi0Average = 0
                if (xi2<=-xi0Trapped) or (xi1>=xi0Trapped) or ((xi1<=-xi0Trapped) and (xi2>=xi0Trapped)):  
                    xi0Average = xi[j]

                if xi0Average != 0:
                    integrand[ir,j,:] =  2*np.pi * p**2 * v * xi0Average / Vprime_VpVol[ir,j,:]

        return integrand
