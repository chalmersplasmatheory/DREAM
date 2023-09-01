# Settings for the runaway electron density

import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
from . PrescribedInitialParameter import PrescribedInitialParameter
from .. import AdvectionInterpolation
from .. TransportSettings import TransportSettings
from . DistributionFunction import DISTRIBUTION_MODE_NUMERICAL



DREICER_RATE_DISABLED = 1
DREICER_RATE_CONNOR_HASTIE_NOCORR= 2
DREICER_RATE_CONNOR_HASTIE = 3
DREICER_RATE_NEURAL_NETWORK = 4

COLLQTY_ECEFF_MODE_EC_TOT = 1
COLLQTY_ECEFF_MODE_CYLINDRICAL = 2
COLLQTY_ECEFF_MODE_SIMPLE = 3
COLLQTY_ECEFF_MODE_FULL = 4

AVALANCHE_MODE_NEGLECT = 1
AVALANCHE_MODE_FLUID = 2
AVALANCHE_MODE_FLUID_HESSLOW = 3
AVALANCHE_MODE_KINETIC = 4

COMPTON_MODE_NEGLECT = 1
COMPTON_MODE_FLUID   = 2
COMPTON_MODE_KINETIC = 3 
COMPTON_RATE_ITER_DMS = -1
COMPTON_RATE_ITER_DMS_KINETIC = -2
ITER_PHOTON_FLUX_DENSITY = 1e18

TRITIUM_MODE_NEGLECT = 1
TRITIUM_MODE_FLUID = 2
TRITIUM_MODE_KINETIC = 3

# Interpolation methods for advection term in transport equation
AD_INTERP_CENTRED  = AdvectionInterpolation.AD_INTERP_CENTRED
AD_INTERP_UPWIND   = AdvectionInterpolation.AD_INTERP_UPWIND
AD_INTERP_UPWIND_2ND_ORDER = AdvectionInterpolation.AD_INTERP_UPWIND_2ND_ORDER
AD_INTERP_DOWNWIND = AdvectionInterpolation.AD_INTERP_DOWNWIND
AD_INTERP_QUICK    = AdvectionInterpolation.AD_INTERP_QUICK 
AD_INTERP_SMART    = AdvectionInterpolation.AD_INTERP_SMART 
AD_INTERP_MUSCL    = AdvectionInterpolation.AD_INTERP_MUSCL 
AD_INTERP_OSPRE    = AdvectionInterpolation.AD_INTERP_OSPRE 
AD_INTERP_TCDF     = AdvectionInterpolation.AD_INTERP_TCDF  

AD_INTERP_JACOBIAN_LINEAR = AdvectionInterpolation.AD_INTERP_JACOBIAN_LINEAR
AD_INTERP_JACOBIAN_FULL   = AdvectionInterpolation.AD_INTERP_JACOBIAN_FULL  
AD_INTERP_JACOBIAN_UPWIND = AdvectionInterpolation.AD_INTERP_JACOBIAN_UPWIND


HOTTAIL_MODE_DISABLED = 1
HOTTAIL_MODE_ANALYTIC = 2 # not yet implemented
HOTTAIL_MODE_ANALYTIC_ALT_PC = 3


class RunawayElectrons(UnknownQuantity,PrescribedInitialParameter):

    def __init__(self, settings, density=0, radius=0, avalanche=AVALANCHE_MODE_NEGLECT, dreicer=DREICER_RATE_DISABLED, compton=COMPTON_MODE_NEGLECT, Eceff=COLLQTY_ECEFF_MODE_FULL, pCutAvalanche=0, comptonPhotonFlux=0, tritium=TRITIUM_MODE_NEGLECT, hottail=HOTTAIL_MODE_DISABLED):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.avalanche = avalanche
        self.dreicer   = dreicer
        self.Eceff     = Eceff
        self.pCutAvalanche = pCutAvalanche
        self.tritium   = tritium
        self.hottail   = hottail
        self.negative_re = False
        self.extrapolateDreicer = True

        self.setCompton(compton, comptonPhotonFlux)

        self.advectionInterpolation = AdvectionInterpolation.AdvectionInterpolation(kinetic=False)
        self.transport = TransportSettings(kinetic=False)

        self.density = None
        self.radius  = None
        self.setInitialProfile(density=density, radius=radius)


    def setInitialProfile(self, density, radius=0):
        _data, _rad = self._setInitialData(data=density, radius=radius)

        self.density = _data
        self.radius  = _rad
        self.verifySettingsPrescribedInitialData()


    def setAvalanche(self, avalanche, pCutAvalanche=0):
        """
        Enables/disables avalanche generation.
        """
        if avalanche == False:
            self.avalanche = AVALANCHE_MODE_NEGLECT
        else:
            self.avalanche = int(avalanche)
            self.pCutAvalanche = pCutAvalanche


    def setDreicer(self, dreicer):
        """
        Specifies which model to use for calculating the
        Dreicer runaway rate.
        """
        if dreicer == False:
            self.dreicer = DREICER_RATE_DISABLED
        else:
            self.dreicer = int(dreicer)


    def setCompton(self, compton, photonFlux = None, photonFlux_t = np.array([0])):
        """
        Specifies which model to use for calculating the
        compton runaway rate.
        """
        if compton == False or compton == COMPTON_MODE_NEGLECT:
            self.compton = COMPTON_MODE_NEGLECT
        else:
            if compton == COMPTON_RATE_ITER_DMS:
                # set fluid compton source and standard ITER flux of 1e18
                compton = COMPTON_MODE_FLUID
                if photonFlux is None:
                    photonFlux = np.array([ITER_PHOTON_FLUX_DENSITY])
                    photonFlux_t = np.array([0])

            if compton == COMPTON_RATE_ITER_DMS_KINETIC:
                # set fluid compton source and standard ITER flux of 1e18
                compton = COMPTON_MODE_KINETIC
                if photonFlux is None:
                    photonFlux = np.array([ITER_PHOTON_FLUX_DENSITY])
                    photonFlux_t = np.array([0])
            
            if photonFlux is None:
                raise EquationException("n_re: Compton photon flux must be set.")
            elif type(photonFlux) == int or type(photonFlux) == float or type(photonFlux) == np.float64:
                photonFlux = np.array([float(photonFlux)])

            self.compton = int(compton)
            self.comptonPhotonFlux = photonFlux
            self.comptonPhotonFlux_t = photonFlux_t


    def setEceff(self, Eceff):
        """
        Specifies which model to use for calculating the
        effective critical field (used in the avalanche formula).
        """
        self.Eceff = int(Eceff)


    def setTritium(self, tritium):
        """
        Specifices whether or not to include runaway generation
        through tritium decay as a source term.
        """
        if tritium == True or tritium == TRITIUM_MODE_FLUID:
            self.tritium = TRITIUM_MODE_FLUID
        if tritium == False or tritium == TRITIUM_MODE_NEGLECT:
            self.tritium = TRITIUM_MODE_NEGLECT
        self.tritium = int(tritium)


    def setHottail(self, hottail):
        """
        Specify which model to use for hottail runaway generation
        """
        if hottail == False:
            self.hottail = HOTTAIL_MODE_DISABLED
        else:
            self.hottail = hottail
            if hottail != HOTTAIL_MODE_DISABLED:
                self.settings.eqsys.f_hot.enableAnalyticalDistribution()


    def setNegativeRunaways(self, negative_re=True):
        """
        Introduce a density of runaway electrons with negative pitch,
        allowing the kinetic avalanche source term to properly account for
        large-angle collisions with runaways moving in different directions.
        """
        self.negative_re = negative_re
        
    def setExtrapolateDreicer(self, extrapolateDreicer=False):
        """
        Extrapolates the result from the neural network for small electric fields
        such that the Dreicer generation rate is continuous and has continuous derivative.
        """
        self.extrapolateDreicer = extrapolateDreicer

    def setAdvectionInterpolationMethod(self, ad_int=AD_INTERP_CENTRED,
        ad_jac=AD_INTERP_JACOBIAN_FULL, fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used in the advection terms of
        the transport equation.
        
        :param int ad_int:               Interpolation method to use for the radial coordinate.
        :param int ad_jac:               Jacobian interpolation mode to use for the radial coordinate.
        :param float fluxlimiterdamping: Damping parameter used to under-relax the interpolation coefficients during non-linear iterations (should be between 0 and 1).
        """
        self.advectionInterpolation.setMethod(ad_int=ad_int, ad_jac=ad_jac, fluxlimiterdamping=fluxlimiterdamping)


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.avalanche = int(data['avalanche'])

        if 'pCutAvalanche' in data:
            self.pCutAvalanche = data['pCutAvalanche']

        self.dreicer   = int(data['dreicer'])
        self.Eceff     = int(data['Eceff'])
        self.compton            = int(data['compton']['mode'])
        self.density   = data['init']['x']
        self.radius    = data['init']['r']

        if 'flux' in data['compton']:
            if type(data['compton']['flux']) == dict:
                self.comptonPhotonFlux  = data['compton']['flux']['x']
                self.comptonPhotonFlux_t = data['compton']['flux']['t']
            else:
                self.comptonPhotonFlux  = data['compton']['flux']
                self.comptonPhotonFlux_t = np.array([0.0])

        if 'adv_interp' in data:
            self.advectionInterpolation.fromdict(data['adv_interp'])

        if 'hottail' in data:
            self.hottail = int(data['hottail'])

        if 'tritium' in data:
            self.tritium = int(data['tritium'])

        if 'negative_re' in data:
            self.negative_re = bool(data['negative_re'])
        
        if 'extrapolateDreicer' in data:
            self.ExtrapolateDreicer = bool(data['extrapolateDreicer'])

        if 'transport' in data:
            self.transport.fromdict(data['transport'])


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this RunawayElectrons object.
        """
        data = {
            'avalanche': self.avalanche,
            'dreicer': self.dreicer,
            'Eceff': self.Eceff,
            'pCutAvalanche': self.pCutAvalanche,
            'transport': self.transport.todict(),
            'tritium': self.tritium,
            'hottail': self.hottail,
            'negative_re': self.negative_re,
            'extrapolateDreicer': self.extrapolateDreicer
        }
        data['compton'] = {
            'mode': self.compton
        }
        if self.compton != COMPTON_MODE_NEGLECT:
            data['compton']['flux'] = {
                'x': self.comptonPhotonFlux,
                't': self.comptonPhotonFlux_t
            }
        data['init'] = {
            'x': self.density,
            'r': self.radius
        }

        # Flux limiter settings
        data['adv_interp'] = self.advectionInterpolation.todict()

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if type(self.avalanche) != int:
            raise EquationException("n_re: Invalid value assigned to 'avalanche'. Expected integer.")
        if type(self.dreicer) != int:
            raise EquationException("n_re: Invalid value assigned to 'dreicer'. Expected integer.")
        if type(self.compton) != int:
            raise EquationException("n_re: Invalid value assigned to 'compton'. Expected integer.")
        if type(self.hottail) != int:
            raise EquationException("n_re: Invalid value assigned to 'hottail'. Expected integer.")
        if type(self.Eceff) != int:
            raise EquationException("n_re: Invalid value assigned to 'Eceff'. Expected integer.")
        if self.avalanche == AVALANCHE_MODE_KINETIC and self.pCutAvalanche == 0:
            raise EquationException("n_re: Invalid value assigned to 'pCutAvalanche'. Must be set explicitly when using KINETIC avalanche.")
        if type(self.tritium) != int:
            raise EquationException("n_re: Invalid value assigned to 'tritium'. Expected integer.")
        if self.hottail != HOTTAIL_MODE_DISABLED and self.settings.eqsys.f_hot.mode == DISTRIBUTION_MODE_NUMERICAL:
            raise EquationException("n_re: Invalid setting combination: when hottail is enabled, the 'mode' of f_hot cannot be NUMERICAL. Enable ANALYTICAL f_hot distribution or disable hottail.")
        if type(self.negative_re) != bool:
            raise EquationException("n_re: Invalid value assigned to 'negative_re'. Expected bool.")
        if type(self.extrapolateDreicer) != bool:
            raise EquationException("n_re: Invalid value assigned to 'extrapolateDreicer'. Expected bool.")

        if self.compton != COMPTON_MODE_NEGLECT:
            if type(self.comptonPhotonFlux) != np.ndarray:
                raise EquationException("Invalid type for 'comptonPhotonFlux'. Expected numpy array.")
            elif type(self.comptonPhotonFlux_t) != np.ndarray:
                raise EquationException("Invalid type for 'comptonPhotonFlux_t'. Expected number array.")
            elif self.comptonPhotonFlux.shape != self.comptonPhotonFlux_t.shape:
                raise EquationException("The shapes of 'comptonPhotonFlux' and 'photonFlux_t' do not match.")

        self.advectionInterpolation.verifySettings()
        self.transport.verifySettings()


    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('n_re', data=self.density, radius=self.radius)


