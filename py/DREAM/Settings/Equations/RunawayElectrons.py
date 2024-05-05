# Settings for the runaway electron density

import numpy as np
from scipy.integrate import quad
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
COMPTON_MACHINE_ITER = 1
ITER_PHOTON_FLUX_DENSITY = 1e18
C1_COMPTON = 1.2
C2_COMPTON = 0.8
C3_COMPTON = 0.

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

# Edge scrape-off term
LCFS_LOSS_MODE_DISABLED = 1
LCFS_LOSS_MODE_FLUID = 2
LCFS_LOSS_MODE_KINETIC = 3


def GammafluxProfil(E, C1, C2, C3):
    mc2 = 0.51099895000e6
    z = (np.log(mc2 * E / 1e6) + C1) / C2 + C3 * (mc2 * E / 1e6)**2
    Gamma = np.exp(-np.exp(-z) - z + 1)
    return Gamma 

class RunawayElectrons(UnknownQuantity,PrescribedInitialParameter):

  
    def __init__(self, settings, density=0, radius=0, avalanche=AVALANCHE_MODE_NEGLECT, dreicer=DREICER_RATE_DISABLED, compton=COMPTON_MODE_NEGLECT, Eceff=COLLQTY_ECEFF_MODE_FULL, pCutAvalanche=0, comptonPhotonFlux=0, C1_Compton=C1_COMPTON, C2_compton=C2_COMPTON, C3_Compton=C3_COMPTON, tritium=TRITIUM_MODE_NEGLECT, hottail=HOTTAIL_MODE_DISABLED, lcfs_loss=LCFS_LOSS_MODE_DISABLED):
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
        
        self.setCompton(compton=compton, photonFlux=comptonPhotonFlux, C1=C1_Compton, C2=C2_compton, C3=C3_Compton)
        self.advectionInterpolation = AdvectionInterpolation.AdvectionInterpolation(kinetic=False)
        self.transport = TransportSettings(kinetic=False)

        self.density = None
        self.radius  = None
        self.setInitialProfile(density=density, radius=radius)
        
        # Loss term
        self.lcfs_loss = lcfs_loss
        self.lcfs_t_loss = 0
        self.lcfs_t_loss_r = 0
        self.lcfs_user_input_psi = 0
        self.lcfs_psi_edge_t0 = 0


    def setInitialProfile(self, density, radius=0):
        _data, _rad = self._setInitialData(data=density, radius=radius)

        self.density = _data
        self.radius  = _rad
        self.verifySettingsPrescribedInitialData()


    # Loss term    
    def setLCFSLoss(self, lcfs_loss):
        """
        Specifies which model to use for calculating the
        LCFS loss term.
        """
        if lcfs_loss == False:
            self.lcfs_loss = LCFS_LOSS_MODE_DISABLED
        else:
            self.lcfs_loss = int(lcfs_loss)
    
    
    # Loss term
    def setLCFSLossTime(self, t_loss, radius=0):
        """
        Sets the timescale constant t_loss for the
        LCFS loss term. 
        """
        _data, _rad = self._setInitialData(data=t_loss, radius=radius)

        self.lcfs_t_loss = _data
        self.lcfs_t_loss_r  = _rad
        self.verifySettingsPrescribedInitialData()
        
        
    def setLCFSLossPsiEdget0(self, psi_edge_t0, user_input_active=1):
        """
        Sets the value of psi_p at the plasma edge at
        t = 0, used to determine the LCFS radial point.
        Use in restarts to keep the value from the first
        simulation. user_input_active should be 1 (default)
        or set to 0 to manually switch off the user input.
        """
        self.lcfs_user_input_psi = int(user_input_active)
        self.lcfs_psi_edge_t0 = psi_edge_t0


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


    def setCompton(self, compton, photonFlux = None, photonFlux_t = np.array([0]), machine = None, C1 = None, C2 = None, C3 = None):
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
            
            if machine == COMPTON_MACHINE_ITER:
                C1 = 1.2
                C2 = 0.8
                C3 = 0.
            else:
                if C1 is None:
                    C1 = C1_COMPTON
                if C2 is None:
                    C2 = C2_COMPTON
                if C3 is None:
                    C3 = C3_COMPTON


            self.compton = int(compton)
            self.comptonPhotonFlux = photonFlux
            self.comptonPhotonFlux_t = photonFlux_t
            self.C1_Compton = C1
            self.C2_Compton = C2
            self.C3_Compton = C3
            self.integratedComptonSpectrum = quad(GammafluxProfil, 0, np.inf, args=(C1, C2, C3))[0]

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
        self.compton   = int(data['compton']['mode'])
        self.density   = data['init']['x']
        self.radius    = data['init']['r']
        # Loss term
        if 'lcfs_loss' in data:
            self.lcfs_loss     = int(data['lcfs_loss'])
            self.lcfs_user_input_psi     = int(data['lcfs_user_input_psi'])
            self.lcfs_psi_edge_t0        = data['lcfs_psi_edge_t0']
            self.lcfs_t_loss   = data['lcfs_t_loss']['x']
            self.lcfs_t_loss_r = data['lcfs_t_loss']['r']

        if 'flux' in data['compton']:
            if type(data['compton']['flux']) == dict:
                self.comptonPhotonFlux  = data['compton']['flux']['x']
                self.comptonPhotonFlux_t = data['compton']['flux']['t']
            else:
                self.comptonPhotonFlux  = data['compton']['flux']
                self.comptonPhotonFlux_t = np.array([0.0])

        if 'C1' in data['compton']:
            self.C1_Compton = float(data['compton']['C1'])
            self.C2_Compton = float(data['compton']['C2'])
            self.C3_Compton = float(data['compton']['C3'])
            self.integratedComptonSpectrum = quad(GammafluxProfil, 0, np.inf, args=(self.C1_Compton, self.C2_Compton, self.C3_Compton))[0]

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
            'lcfs_loss': self.lcfs_loss, # Loss term
            'lcfs_user_input_psi': self.lcfs_user_input_psi,
            'lcfs_psi_edge_t0': self.lcfs_psi_edge_t0,
            'extrapolateDreicer': self.extrapolateDreicer
        }
        if self.compton != COMPTON_MODE_NEGLECT:
            data['compton'] = {
                'mode': self.compton, 
                'C1': self.C1_Compton,
                'C2': self.C2_Compton,
                'C3': self.C3_Compton,
                'gammaInt': self.integratedComptonSpectrum
            }
            data['compton']['flux'] = {
                'x': self.comptonPhotonFlux,
                't': self.comptonPhotonFlux_t
            }
        else:
            data['compton'] = {
                'mode': self.compton
            }
        data['init'] = {
            'x': self.density,
            'r': self.radius
        }
        # Loss term
        data['lcfs_t_loss'] = {
            'x': self.lcfs_t_loss,
            'r': self.lcfs_t_loss_r
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
            if type(self.integratedComptonSpectrum) != float:
                raise EquationException("Invalid type for 'integratedComptonSpectrum'. Expected float.")
            if type(self.C1_Compton) != float:
                raise EquationException("Invalid type for 'C1_Compton'. Expected float.")
            if type(self.C2_Compton) != float:
                raise EquationException("Invalid type for 'C2_Compton'. Expected float.")
            if type(self.C3_Compton) != float:
                raise EquationException("Invalid type for 'C3_Compton'. Expected float.")

        self.advectionInterpolation.verifySettings()
        self.transport.verifySettings()


    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('n_re', data=self.density, radius=self.radius)


