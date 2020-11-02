
import numpy as np

from .. Settings.MomentumGrid import TYPE_PXI, TYPE_PPARPPERP


class MomentumGrid:
    

    def __init__(self, name, rgrid, data):
        """
        Constructor.

        name:  Grid name.
        rgrid: Parent 'Grid' object (representing radial grid).
        data:  Momentum grid data.
        """
        self.name  = name
        self.r     = rgrid.r
        self.r_f   = rgrid.r_f
        self.dr    = rgrid.dr
        self.type  = data['type']
        self.rgrid = rgrid

        self.Vprime = data['Vprime']
        self.p1   = data['p1']
        self.p2   = data['p2']
        self.p1_f = data['p1_f']
        self.p2_f = data['p2_f']
        self.dp1  = data['dp1']
        self.dp2  = data['dp2']

        self.Vprime_VpVol = np.copy(self.Vprime)
        for i in range(0, self.rgrid.r.size):
            self.Vprime_VpVol[i,:] /= rgrid.VpVol[i]

        self.DR, self.DP2, self.DP1 = np.meshgrid(self.dr, self.dp2, self.dp1, indexing='ij')
        #self.DP2, self.DR, self.DP1 = np.meshgrid(self.dp2, self.dr, self.dp1)


    def integrate3D(self, data, axes=(-3,-2,-1)):
        """
        Evaluate a numerical volume integral of the given data
        on this grid.
        
        data: Data to numerically integrate.
        axes: Axes to integrate over.
        """
        if len(axes) != 3:
            raise OutputException("Invalid 'axes' parameter provided to 'integrate3D()'.")

        return (data * self.Vprime * self.DR * self.DP1 * self.DP2).sum(axes)


    def integrate2D(self, data, axes=(-2,-1)):
        """
        Evaluate a numerical momentum integral of the given
        data on this grid.

        data: Data to numerically integrate.
        axes: Axes to integrate over.
        """
        if len(axes) != 2:
            raise OutputException("Invalid 'axes' parameter provided to 'integrate2D()'.")

        return (data * (self.Vprime/self.rgrid.VpVol[:,None,None]) * self.DP1 * self.DP2).sum(axes)


    def getGamma(self):
        """
        Returns a meshgrid representing the relativistic factor on
        this 2D momentum grid.
        """
        raise OutputException("'getGamma()' has not been implemented for this momentum grid.")


    def getP1TeXName(self):
        """
        Returns the TeX-compatible name of the p1 coordinate.
        """
        if self.type == TYPE_PXI:
            return r'$p$'
        elif self.type == TYPE_PPARPPERP:
            return r'$p_\parallel$'
        else:
            raise OutputException("Unrecognized grid type: {}".format(self.type))
            

    def getP2TeXName(self):
        """
        Returns the TeX-compatible name of the p2 coordinate.
        """
        if self.type == TYPE_PXI:
            return r'$\xi$'
        elif self.type == TYPE_PPARPPERP:
            return r'$p_\perp$'
        else:
            raise OutputException("Unrecognized grid type: {}".format(self.type))


    def getVpar(self):
        """
        Returns a meshgrid representing the parallel velocity on this
        2D momentum grid. (This method must be implemented separately
        for each specific momentum grid type)
        """
        raise OutputException("'getVpar()' has not been implemented for this momentum grid.")


