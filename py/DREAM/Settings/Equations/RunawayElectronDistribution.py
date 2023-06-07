
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from . import DistributionFunction as DistFunc
from . DistributionFunction import DistributionFunction
from .. TransportSettings import TransportSettings


INIT_FORWARD = 1
INIT_XI_NEGATIVE = 2
INIT_XI_POSITIVE = 3
INIT_ISOTROPIC = 4
INIT_AVALANCHE = 5


class RunawayElectronDistribution(DistributionFunction):
    
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


    def setInitType(self, inittype):
        """
        Specifies how the runaway electron distribution function f_re should be
        initialized from the runaway density n_re.

        :param int inittype: Flag indicating how to initialize f_re.
        """
        self.inittype = int(inittype)


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


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this RunawayElectronDistribution object.
        """
        d = super().todict()
        d['inittype'] = self.inittype

        return d


