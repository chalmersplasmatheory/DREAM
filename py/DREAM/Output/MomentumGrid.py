
import numpy as np


class MomentumGrid:
    

    def __init__(self, name, r, dr, data):
        """
        Constructor.

        name: Grid name.
        r:    Radial grid.
        data: Momentum grid data.
        """
        self.name = name
        self.r    = r
        self.dr   = dr
        self.type = data['type']

        self.Vprime = data['Vprime']
        self.p1  = data['p1']
        self.p2  = data['p2']
        self.dp1 = data['dp1']
        self.dp2 = data['dp2']

        #self.DR, self.DP1, self.DP2 = np.meshgrid(self.r, self.p1, self.p2)


    def trapz(self, data):
        """
        Evaluate a numerical volume integral of the given data
        using a trapezoidal rule on this grid.
        """
        #I = data * self.Vprime * self.DR * self.DP1 * self.DP2
        I = np.trapz(data * Vprime, dx=self.dp2, axis=2)
        I = np.trapz(I, dx=self.dp1, axis=1)
        I = np.trapz(I, dx=self.dr, axis=0)

        return I


