
import numpy as np
import scipy.constants

from .OutputException import OutputException

class MomentumGrid:
    P1_NAME = None
    P2_NAME = None
    P1_TEX_NAME = None
    P2_TEX_NAME = None
    
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
        if 'Vprime_f2' in data:
            self.Vprime_f2 = data['Vprime_f2']
        self.p1   = data['p1']
        self.p2   = data['p2']
        self.p1_f = data['p1_f']
        self.p2_f = data['p2_f']
        self.dp1  = data['dp1']
        self.dp2  = data['dp2']

        self.Vprime_VpVol = np.copy(self.Vprime[:])
        for i in range(0, self.rgrid.r.size):
            self.Vprime_VpVol[i,:] /= rgrid.VpVol[i]

        self.DR, self.DP2, self.DP1 = np.meshgrid(self.dr[:], self.dp2[:], self.dp1[:], indexing='ij')

        # populated by subclass constructors
        self._PPAR  = None
        self._PPERP = None
        self._P     = None
        self._XI    = None
        self._GAMMA = None

        self._PPAR_f  = None
        self._PPERP_f = None
        self._P_f     = None
        self._XI_f    = None

        self._VprimeCylindrical = None
        self._VprimeSpherical   = None

    @property
    def p1name(self):
        if self.P1_NAME is None:
            raise OutputException("P1_NAME has not been implemented.")
        return self.P1_NAME

    @property
    def p2name(self):
        if self.P2_NAME is None:
            raise OutputException("P2_NAME has not been implemented.")
        return self.P2_NAME

    @property
    def p1TeXname(self):
        if self.P1_TEX_NAME is None:
            raise OutputException("P1_TEX_NAME has not been implemented.")
        return self.P1_TEX_NAME

    @property
    def p2TeXname(self):
        if self.P2_TEX_NAME is None:
            raise OutputException("P2_TEX_NAME has not been implemented.")
        return self.P2_TEX_NAME

    @property
    def P(self):
        if self._P is None:
            raise OutputException("P calculation has not been implemented.")
        return self._P

    @property
    def XI(self):
        if self._XI is None:
            raise OutputException("XI calculation has not been implemented.")
        return self._XI

    @property
    def PPAR(self):
        if self._PPAR is None:
            raise OutputException("PPAR calculation has not been implemented.")
        return self._PPAR

    @property
    def PPERP(self):
        if self._PPERP is None:
            raise OutputException("PPERP calculation has not been implemented.")
        return self._PPERP

    @property
    def P_f(self):
        if self._P_f is None:
            raise OutputException("P_f calculation has not been implemented.")
        return self._P_f

    @property
    def XI_f(self):
        if self._XI_f is None:
            raise OutputException("XI_f calculation has not been implemented.")
        return self._XI_f

    @property
    def PPAR_f(self):
        if self._PPAR_f is None:
            raise OutputException("PPAR_f calculation has not been implemented.")
        return self._PPAR_f

    @property
    def PPERP_f(self):
        if self._PPERP_f is None:
            raise OutputException("PPERP_f calculation has not been implemented.")
        return self._PPERP_f

    @property
    def GAMMA(self):
        if self._GAMMA is None:
            raise OutputException("GAMMA calculation has not been implemented.")
        return self._GAMMA

    @property
    def VprimeCylindrical(self):
        if self._VprimeCylindrical is None:
            raise OutputException("VprimeCylindrical calculation has not been implemented.")
        return self._VprimeCylindrical

    @property
    def VprimeSpherical(self):
        if self._VprimeSpherical is None:
            raise OutputException("VprimeSpherical calculation has not been implemented.")
        return self._VprimeSpherical


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
        return (data * self.Vprime_VpVol * self.DP1 * self.DP2).sum(axes)


    def getGamma(self):
        """
        Returns a meshgrid representing the relativistic factor on
        this 2D momentum grid.
        """
        return self.GAMMA


    def getP1TeXName(self):
        """
        Returns the TeX-compatible name of the p1 coordinate.
        """
        return self.p1TeXname            


    def getP2TeXName(self):
        """
        Returns the TeX-compatible name of the p2 coordinate.
        """
        return self.p2TeXname


    def getVpar(self):
        """
        Returns a meshgrid representing the parallel velocity on this
        2D momentum grid. (This method must be implemented separately
        for each specific momentum grid type)
        """
        return scipy.constants.c * (self.PPAR/self.GAMMA)


    def getBounceAveragedVpar(self):
        """
        Returns a meshgrid representing the integrand that should weigh
        a function when carrying out the v_par moment of a quantity.

        It should be identical to the ``integrand`` produced by the 
        ``CurrentDensityFromDistributionFunction`` class in the DREAM kernel.
        
        (This method must be implemented separately for each specific momentum grid type)
        """
        raise OutputException("'getBounceAveragedVpar()' has not been implemented for this momentum grid.")


