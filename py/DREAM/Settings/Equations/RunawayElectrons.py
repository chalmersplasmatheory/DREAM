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

COLLQTY_ECEFF_MODE_CYLINDRICAL = 1
COLLQTY_ECEFF_MODE_SIMPLE = 2
COLLQTY_ECEFF_MODE_FULL = 3

AVALANCHE_MODE_NEGLECT = 1
AVALANCHE_MODE_FLUID = 2
AVALANCHE_MODE_FLUID_HESSLOW = 3
AVALANCHE_MODE_KINETIC = 4

COMPTON_MODE_NEGLECT = 1
COMPTON_MODE_FLUID   = 2
COMPTON_MODE_KINETIC = 3 
COMPTON_RATE_ITER_DMS = -1
ITER_PHOTON_FLUX_DENSITY = 1e18

class RunawayElectrons(UnknownQuantity,PrescribedInitialParameter):
    

    def __init__(self, settings, density=0, radius=0, avalanche=AVALANCHE_MODE_NEGLECT, dreicer=DREICER_RATE_DISABLED, compton=COMPTON_MODE_NEGLECT, Eceff=COLLQTY_ECEFF_MODE_CYLINDRICAL, pCutAvalanche=0, comptonPhotonFlux=0):
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


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        #self.avalanche = data['avalanche']
        self.avalanche = int(data['avalanche'])
        self.pCutAvalanche = data['pCutAvalanche']
        #self.dreicer   = data['dreicer']
        self.dreicer   = int(data['dreicer'])
        self.Eceff     = int(data['Eceff'])
        self.compton            = int(data['compton']['mode'])
        self.comptonPhotonFlux  = data['compton']['flux']
        self.density   = data['init']['x']
        self.radius    = data['init']['r']

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
            'transport': self.transport.todict()
        }
        data['compton'] = {
            'mode': self.compton,
            'flux': self.comptonPhotonFlux
        }
        data['init'] = {
            'x': self.density,
            'r': self.radius
        }

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



    def verifySettingsPrescribedInitialData(self):
        self._verifySettingsPrescribedInitialData('n_re', data=self.density, radius=self.radius)
