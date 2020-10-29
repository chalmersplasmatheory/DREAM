
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from . import DistributionFunction as DistFunc
from . DistributionFunction import DistributionFunction
from .. TransportSettings import TransportSettings


class RunawayElectronDistribution(DistributionFunction):
    
    def __init__(self, settings,
        fre=None, initr=None, initp=None, initxi=None,
        initppar=None, initpperp=None,
        rn0=None, n0=None, rT0=None, T0=None, bc=DistFunc.BC_PHI_CONST,
        ad_int_r=DistFunc.AD_INTERP_CENTRED,
        ad_int_p1=DistFunc.AD_INTERP_CENTRED,
        ad_int_p2=DistFunc.AD_INTERP_CENTRED,
        fluxlimiterdamping=1.0):
        """
        Constructor.
        """
        super().__init__(settings=settings, name='f_re', grid=settings.runawaygrid,
            f=fre, initr=initr, initp=initp, initxi=initxi, initppar=initppar,
            initpperp=initpperp, rn0=rn0, n0=n0, rT0=rT0, T0=T0,
            bc=bc, ad_int_r=ad_int_r, ad_int_p1=ad_int_p1,
            ad_int_p2=ad_int_p2, fluxlimiterdamping=fluxlimiterdamping)


    def fromdict(self, data):
        """
        Load data for this object from the given dictionary.
        """
        super().fromdict(data)


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this RunawayElectronDistribution object.
        """
        return super().todict()


