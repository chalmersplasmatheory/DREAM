
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from . import DistributionFunction as DistFunc
from . DistributionFunction import DistributionFunction
from . PrescribedInitialParameter import PrescribedInitialParameter
from .. TransportSettings import TransportSettings


INIT_FORWARD = 1
INIT_XI_NEGATIVE = 2
INIT_XI_POSITIVE = 3
INIT_ISOTROPIC = 4
INIT_AVALANCHE = 5
INIT_PRESCRIBED = 6


class RunawayElectronDistribution(DistributionFunction, PrescribedInitialParameter):
    
    def __init__(self, settings,
        fre=[0.0], initr=[0.0], initp=[0.0], initxi=[0.0],
        initppar=None, initpperp=None,
        rn0=None, n0=None, rT0=None, T0=None, bc=DistFunc.BC_PHI_CONST,
        ad_int_r=DistFunc.AD_INTERP_CENTRED,
        ad_int_p1=DistFunc.AD_INTERP_CENTRED,
        ad_int_p2=DistFunc.AD_INTERP_CENTRED,
        ad_jac_r=DistFunc.AD_INTERP_JACOBIAN_LINEAR,
        ad_jac_p1=DistFunc.AD_INTERP_JACOBIAN_LINEAR, 
        ad_jac_p2=DistFunc.AD_INTERP_JACOBIAN_LINEAR,
        fluxlimiterdamping=1.0):
        """
        Constructor.
        """
        super().__init__(settings=settings, name='f_re', grid=settings.runawaygrid,
            f=fre, initr=initr, initp=initp, initxi=initxi, initppar=initppar,
            initpperp=initpperp, rn0=rn0, n0=n0, rT0=rT0, T0=T0,
            bc=bc, ad_int_r=ad_int_r, ad_int_p1=ad_int_p1,
            ad_int_p2=ad_int_p2, fluxlimiterdamping=fluxlimiterdamping)

        self.inittype = INIT_FORWARD

    def setInitialValue(self, f, r, p=None, xi=None, ppar=None, pperp=None):
        """
        Set the initial value of this electron distribution function. Only one
        of the pairs (p, xi) and (ppar, pperp) of momentum grids need to be
        given.

        :param f:     Array representing the distribution function value on the grid (must have size (nr, nxi, np) or (nr, npperp, nppar))
        :param r:     Radial grid on which the initial distribution is given.
        :param p:     Momentum grid.
        :param xi:    Pitch grid.
        :param ppar:  Parallel momentum grid.
        :param pperp: Perpendicular momentum grid.
        """
        super().setInitialValue(f, r, p=p, xi=xi, ppar=ppar, pperp=pperp)

        self.E_init = 1

        self.setInitType(INIT_PRESCRIBED)


    def setInitType(self, inittype):
        """
        Specifies how the runaway electron distribution function f_re should be
        initialized from the runaway density n_re.

        :param int inittype: Flag indicating how to initialize f_re.
        """
        self.inittype = int(inittype)


    def setInitialAvalancheDistribution(self, E=1, r=0):
        """
        Initialize the runaway distribution function according to an analytical
        avalanche distribution.
        """
        self.setInitType(INIT_AVALANCHE)

        _data, _rad = self._setInitialData(data=E, radius=r)
        self.E_init = _data
        self.E_init_r = _rad

        self.verifySettingsPrescribedInitialData()


    def fromdict(self, data):
        """
        Load data for this object from the given dictionary.
        """
        super().fromdict(data)

        def scal(v):
            if type(v) == np.ndarray: return v[0]
            else: return v

        if 'inittype' in data:
            self.inittype = int(scal(data['inittype']))

        if 'E_init' in data:
            self.E_init = data['E_init']['x']
            self.E_init_r = data['E_init']['r']


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this RunawayElectronDistribution object.
        """
        d = super().todict()
        d['inittype'] = self.inittype

        if self.inittype == INIT_AVALANCHE:
            d['E_init'] = {
                'x': self.E_init,
                'r': self.E_init_r
            }

        return d


    def verifySettings(self):
        """
        Verify the settings of this module.
        """
        self.verifySettingsPrescribedInitialData()


    def verifySettingsPrescribedInitialData(self):
        """
        Verify that the prescribed initial data is correctly set.
        """
        if self.inittype == INIT_AVALANCHE:
            self._verifySettingsPrescribedInitialData('E_init', data=self.E_init, radius=self.E_init_r)


