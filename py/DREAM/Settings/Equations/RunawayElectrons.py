# Settings for the runaway electron density

import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
from . PrescribedInitialParameter import PrescribedInitialParameter
from .. TransportSettings import TransportSettings



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
ITER_PHOTON_FLUX_DENSITY = 1e18

# Interpolation methods for advection term in transport equation
AD_INTERP_CENTRED  = 1
AD_INTERP_UPWIND   = 2
AD_INTERP_UPWIND_2ND_ORDER = 3
AD_INTERP_DOWNWIND = 4
AD_INTERP_QUICK    = 5
AD_INTERP_SMART    = 6
AD_INTERP_MUSCL    = 7
AD_INTERP_OSPRE    = 8
AD_INTERP_TCDF     = 9

AD_INTERP_JACOBIAN_LINEAR = 1
AD_INTERP_JACOBIAN_FULL   = 2
AD_INTERP_JACOBIAN_UPWIND = 3


class RunawayElectrons(UnknownQuantity,PrescribedInitialParameter):
    

    def __init__(self, settings, density=0, radius=0, avalanche=AVALANCHE_MODE_NEGLECT, dreicer=DREICER_RATE_DISABLED, compton=COMPTON_MODE_NEGLECT, Eceff=COLLQTY_ECEFF_MODE_CYLINDRICAL, pCutAvalanche=0, comptonPhotonFlux=0, tritium=False):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.avalanche = avalanche
        self.dreicer   = dreicer
        self.compton   = compton
        self.comptonPhotonFlux = comptonPhotonFlux
        self.Eceff     = Eceff
        self.pCutAvalanche = pCutAvalanche
        self.tritium   = tritium

        self.transport = TransportSettings(kinetic=False)

        self.adv_interp_r  = AD_INTERP_CENTRED 
        self.adv_jac_r  = AD_INTERP_JACOBIAN_FULL
        self.fluxlimiterdamping = 1.0

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
        self.avalanche = int(avalanche)
        self.pCutAvalanche = pCutAvalanche


    def setDreicer(self, dreicer):
        """
        Specifies which model to use for calculating the
        Dreicer runaway rate.
        """
        self.dreicer = int(dreicer)

    def setCompton(self, compton, photonFlux = None):
        """
        Specifies which model to use for calculating the
        compton runaway rate.
        """
        if compton == COMPTON_RATE_ITER_DMS:
            # set fluid compton source and standard ITER flux of 1e18
            compton = COMPTON_MODE_FLUID
            if photonFlux is None:
                photonFlux = ITER_PHOTON_FLUX_DENSITY
        
        if photonFlux is None:
            raise EquationException("n_re: Compton photon flux must be set.")

        self.compton = int(compton)
        self.comptonPhotonFlux = photonFlux

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
        self.tritium = tritium


    def setAdvectionInterpolationMethod(self, ad_int=AD_INTERP_CENTRED,
        ad_jac=AD_INTERP_JACOBIAN_FULL, fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used in the advection terms of
        the kinetic equation. To set all three components, provide ad_int and/or ad_jac.
        Otherwise the three components can use separate interpolation methods.
        
        :param int ad_int:               Interpolation method to use for the radial coordinate.
        :param int ad_jac:               Jacobian interpolation mode to use for the radial coordinate.
        :param float fluxlimiterdamping: Damping parameter used to under-relax the interpolation coefficients during non-linear iterations (should be between 0 and 1).
        """
        self.fluxlimiterdamping = fluxlimiterdamping
        self.adv_interp_r = ad_int
        self.adv_jac_r = ad_jac


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        self.avalanche = int(data['avalanche'])
        self.pCutAvalanche = data['pCutAvalanche']
        self.dreicer   = int(data['dreicer'])
        self.Eceff     = int(data['Eceff'])
        self.compton            = int(data['compton']['mode'])
        self.comptonPhotonFlux  = data['compton']['flux']
        self.density   = data['init']['x']
        self.radius    = data['init']['r']

        if 'adv_interp' in data:
            self.adv_interp_r = data['adv_interp']['r']
            self.fluxlimiterdamping = data['adv_interp']['fluxlimiterdamping']
        if 'adv_jac_mode' in data:
            self.adv_jac_r = data['adv_jac_mode']['r']

        if 'tritium' in data:
            self.tritium = bool(data['tritium'])

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
            'tritium': self.tritium
        }
        data['compton'] = {
            'mode': self.compton,
            'flux': self.comptonPhotonFlux
        }
        data['init'] = {
            'x': self.density,
            'r': self.radius
        }

        # Flux limiter settings
        data['adv_interp'] = {}
        data['adv_interp']['r']  = self.adv_interp_r
        data['adv_jac_mode'] = {}
        data['adv_jac_mode']['r'] = self.adv_jac_r
        data['adv_interp']['fluxlimiterdamping'] = self.fluxlimiterdamping

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
        if type(self.Eceff) != int:
            raise EquationException("n_re: Invalid value assigned to 'Eceff'. Expected integer.")
        if self.avalanche == AVALANCHE_MODE_KINETIC and self.pCutAvalanche == 0:
            raise EquationException("n_re: Invalid value assigned to 'pCutAvalanche'. Must be set explicitly when using KINETIC avalanche.")
        if type(self.tritium) != bool:
            raise EquationException("n_re: Invalid value assigned to 'tritium'. Expected bool.")

        ad_int_opts = [
            AD_INTERP_CENTRED, AD_INTERP_DOWNWIND, AD_INTERP_UPWIND, AD_INTERP_UPWIND_2ND_ORDER, 
            AD_INTERP_QUICK, AD_INTERP_SMART, AD_INTERP_MUSCL, AD_INTERP_OSPRE, AD_INTERP_TCDF
        ]
        if self.adv_interp_r not in ad_int_opts:
            raise EquationException("{}: Invalid radial interpolation coefficient set: {}.".format(self.name, self.adv_interp_r))

        self.transport.verifySettings()


    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('n_re', data=self.density, radius=self.radius)


