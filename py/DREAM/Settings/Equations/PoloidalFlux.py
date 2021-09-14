import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings

class PoloidalFlux(UnknownQuantity,PrescribedParameter):
    
    def __init__(self, settings):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.hyperresistivity_enabled = False
        self.hyperresistivity_Lambda_x = None
        self.hyperresistivity_Lambda_r = None
        self.hyperresistivity_Lambda_t = None


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        if 'hyperresistivity' in data:
            hyp = data['hyperresistivity']
            self.hyperresistivity_enabled = bool(hyp['enabled'])

            if 'Lambda' in hyp:
                self.hyperresistivity_Lambda_x = hyp['Lambda']['x']
                self.hyperresistivity_Lambda_r = hyp['Lambda']['r']
                self.hyperresistivity_Lambda_t = hyp['Lambda']['t']


    def setHyperresistivity(self, Lambda, radius=None, times=None):
        """
        Enable the hyperresistive diffusion term and specify the
        transport coefficient ``Lambda``.

        :param Lambda: Diffusion coefficient.
        :param radius: Radial grid on which  ``Lambda`` is specified (if any).
        :param times:  Time grid on which ``Lambda`` is specified (if any).
        """
        d, r, t = self._setPrescribedData(data=Lambda, radius=radius, times=times)

        self.hyperresistivity_enabled = True
        self.hyperresistivity_Lambda_x = d
        self.hyperresistivity_Lambda_r = r
        self.hyperresistivity_Lambda_t = t


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this PoloidalFlux object.
        """
        data = {
            'hyperresistivity': {
                'enabled': self.hyperresistivity_enabled
            }
        }

        if self.hyperresistivity_enabled:
            data['hyperresistivity']['Lambda'] = {
                'x': self.hyperresistivity_Lambda_x,
                'r': self.hyperresistivity_Lambda_r,
                't': self.hyperresistivity_Lambda_t
            }

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.hyperresistivity_enabled:
            self._verifySettingsPrescribedData('psi_p hyperresistivity', data=self.hyperresistivity_Lambda_x, radius=self.hyperresistivity_Lambda_r, times=self.hyperresistivity_Lambda_t)


