
import numpy as np

from .. Settings.MomentumGrid import MOMENTUMGRID_TYPE_PXI, MOMENTUMGRID_TYPE_PPARPPERP


class MomentumGrid:
    

    def __init__(self, name, r, r_f, dr, data):
        """
        Constructor.

        name: Grid name.
        r:    Radial grid.
        data: Momentum grid data.
        """
        self.name = name
        self.r    = r
        self.r_f  = r_f
        self.dr   = dr
        self.type = data['type']

        self.Vprime = data['Vprime']
        self.p1   = data['p1']
        self.p2   = data['p2']
        self.p1_f = data['p1_f']
        self.p2_f = data['p2_f']
        self.dp1  = data['dp1']
        self.dp2  = data['dp2']

        self.DR, self.DP1, self.DP2 = np.meshgrid(self.dr, self.dp1, self.dp2)


    def integrate(self, data, axes=(-3,-2,-1)):
        """
        Evaluate a numerical volume integral of the given data
        using a trapezoidal rule on this grid.
        
        axes: Axes to integrate over.
        """
        if len(axes) != 3:
            raise OutputException("Invalid 'axes' parameter provided to 'integrate()'.")

        return (data * self.Vprime * self.DR * self.DP1 * self.DP2).sum(axes)


    def getP1TeXName(self):
        """
        Returns the TeX-compatible name of the p1 coordinate.
        """
        if self.type == MOMENTUMGRID_TYPE_PXI:
            return r'$p$'
        elif self.type == MOMENTUMGRID_TYPE_PPARPPERP:
            return r'$p_\parallel$'
        else:
            raise OutputException("Unrecognized grid type: {}".format(self.type))
            

    def getP2TeXName(self):
        """
        Returns the TeX-compatible name of the p2 coordinate.
        """
        if self.type == MOMENTUMGRID_TYPE_PXI:
            return r'$\xi$'
        elif self.type == MOMENTUMGRID_TYPE_PPARPPERP:
            return r'$p_\perp$'
        else:
            raise OutputException("Unrecognized grid type: {}".format(self.type))


