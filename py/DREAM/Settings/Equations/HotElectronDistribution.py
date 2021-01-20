
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from . import DistributionFunction as DistFunc
from . DistributionFunction import DistributionFunction
from .. TransportSettings import TransportSettings


# BOUNDARY CONDITIONS (WHEN f_re IS DISABLED)
# (NOTE: These are kept for backwards compatibility. You
#  should _really_ use 'DistributionFunction.XXX' instead)
BC_F_0        = DistFunc.BC_F_0
BC_PHI_CONST  = DistFunc.BC_PHI_CONST
BC_DPHI_CONST = DistFunc.BC_DPHI_CONST

# Interpolation methods for advection term in kinetic equation
AD_INTERP_CENTRED  = DistFunc.AD_INTERP_CENTRED
AD_INTERP_UPWIND   = DistFunc.AD_INTERP_UPWIND
AD_INTERP_UPWIND_2ND_ORDER = DistFunc.AD_INTERP_UPWIND_2ND_ORDER
AD_INTERP_DOWNWIND = DistFunc.AD_INTERP_DOWNWIND
AD_INTERP_QUICK    = DistFunc.AD_INTERP_QUICK
AD_INTERP_SMART    = DistFunc.AD_INTERP_SMART
AD_INTERP_MUSCL    = DistFunc.AD_INTERP_MUSCL
AD_INTERP_OSPRE    = DistFunc.AD_INTERP_OSPRE
AD_INTERP_TCDF     = DistFunc.AD_INTERP_TCDF

AD_INTERP_JACOBIAN_LINEAR = DistFunc.AD_INTERP_JACOBIAN_LINEAR
AD_INTERP_JACOBIAN_FULL   = DistFunc.AD_INTERP_JACOBIAN_FULL
AD_INTERP_JACOBIAN_UPWIND = DistFunc.AD_INTERP_JACOBIAN_UPWIND

HOT_REGION_P_MODE_MC = 1
HOT_REGION_P_MODE_THERMAL = 2
HOT_REGION_P_MODE_THERMAL_SMOOTH = 3

PARTICLE_SOURCE_ZERO     = 1
PARTICLE_SOURCE_IMPLICIT = 2
PARTICLE_SOURCE_EXPLICIT = 3

class HotElectronDistribution(DistributionFunction):
    
    def __init__(self, settings,
        fhot=None, initr=None, initp=None, initxi=None,
        initppar=None, initpperp=None,
        rn0=None, n0=None, rT0=None, T0=None, bc=BC_PHI_CONST,
        ad_int_r =AD_INTERP_CENTRED,
        ad_int_p1=AD_INTERP_CENTRED,
        ad_int_p2=AD_INTERP_CENTRED,
        ad_jac_r =AD_INTERP_JACOBIAN_FULL,
        ad_jac_p1=AD_INTERP_JACOBIAN_FULL, 
        ad_jac_p2=AD_INTERP_JACOBIAN_FULL,
        fluxlimiterdamping=1.0,
        pThreshold=10, pThresholdMode=HOT_REGION_P_MODE_THERMAL,
        particleSource=PARTICLE_SOURCE_EXPLICIT):
        """
        Constructor.
        """
        super().__init__(settings=settings, name='f_hot', grid=settings.hottailgrid,
            f=fhot, initr=initr, initp=initp, initxi=initxi, initppar=initppar,
            initpperp=initpperp, rn0=rn0, n0=n0, rT0=rT0, T0=T0,
            bc=bc, ad_int_r=ad_int_r, ad_int_p1=ad_int_p1,
            ad_int_p2=ad_int_p2, ad_jac_r=ad_jac_r, ad_jac_p1=ad_jac_p1,
            ad_jac_p2=ad_jac_p2, fluxlimiterdamping=fluxlimiterdamping)

        self.pThreshold     = pThreshold
        self.pThresholdMode = pThresholdMode

        self.particleSource = particleSource


    def setHotRegionThreshold(self, pThreshold=10, pMode=HOT_REGION_P_MODE_THERMAL):
        """
        Sets the boundary 'pThreshold' which defines the cutoff separating 'cold'
        from 'hot' electrons when using collfreq_mode FULL. 
        """
        self.pThreshold = pThreshold
        self.pThresholdMode = pMode

    def setParticleSource(self,particleSource=PARTICLE_SOURCE_EXPLICIT):
        """
        Sets which model to use for S_particle if using collfreq_mode FULL,
        which is designed to force the density moment of f_hot to n_cold+n_hot.

        ZERO: The particle source is disabled and set to zero
        EXPLICIT/IMPLICIT: Two in principle equivalent models, but can be more or less stable in different situations. 
        """
        self.particleSource=particleSource


    def fromdict(self, data):
        """
        Load data for this object from the given dictionary.
        """
        super().fromdict(data)
        if 'pThreshold' in data:
            self.pThreshold = data['pThreshold']
            self.pThresholdMode = data['pThresholdMode']
        if 'particleSource' in data:
            self.particleSource = data['particleSource']


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this HotElectronDistribution object.
        """
        data = super().todict()
        if self.grid.enabled:
            data['pThreshold']     = self.pThreshold
            data['pThresholdMode'] = self.pThresholdMode
            data['particleSource'] = self.particleSource

        return data


