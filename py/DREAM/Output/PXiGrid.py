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
        self.xi_f = data['p2_f']
        self.dp = data['dp1']
        self.dxi = data['dp2']

        self.P, self.XI = np.meshgrid(self.p[:], self.xi[:])
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

    def getBounceAveragedVpar(self):
        """
        Returns a meshgrid representing the integrand that should weigh
        a function when carrying out the v_par moment of a quantity.
        
        It should be identical to the ``integrand`` produced by the 
        ``CurrentDensityFromDistributionFunction`` class in the DREAM kernel.        
        """
        c  = scipy.constants.speed_of_light
        integrand = np.zeros(self.Vprime.shape)
        for ir in range(0, self.rgrid.r.size):
            xi0Trapped = self.rgrid.xi0TrappedBoundary[ir]
            for j in range(0, self.xi.size):
                xi1 = self.xi_f[j]
                xi2 = self.xi_f[j+1]
                if(xi1>xi2):
                    xi_t = xi1
                    xi1 = xi2
                    xi2 = xi_t
                
                xi0Average = 0
                if (xi2<=-xi0Trapped) or (xi1>=xi0Trapped) or ((xi1<=-xi0Trapped) and (xi2>=xi0Trapped)):  
                    xi0Average = self.xi[j]

                if xi0Average != 0:
                    for i in range(0, self.p.size):
                        v = c * self.p[i] / np.sqrt(1+self.p[i]**2)
                        integrand[ir,j,i] =  2*np.pi * self.p[i]**2 * v * xi0Average / self.Vprime_VpVol[ir,j,i]

        return integrand





