import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings


HYPERRESISTIVITY_TYPE_NEGLECT = 1
HYPERRESISTIVITY_TYPE_PRESCRIBED = 2
HYPERRESISTIVITY_TYPE_ADAPTIVE = 3


class PoloidalFlux(UnknownQuantity,PrescribedParameter):
    
    def __init__(self, settings):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.hyperresistivity_type = HYPERRESISTIVITY_TYPE_PRESCRIBED
        self.hyperresistivity_Lambda_x = None
        self.hyperresistivity_Lambda_r = None
        self.hyperresistivity_Lambda_t = None

        self.hyperresistivity_grad_j_tot_max = None
        self.hyperresistivity_Lambda0 = None
        self.hyperresistivity_min_duration = None


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        if 'hyperresistivity' in data:
            hyp = data['hyperresistivity']
            self.hyperresistivity_type = int(hyp['type'])

            if 'Lambda' in hyp:
                self.hyperresistivity_Lambda_x = hyp['Lambda']['x']
                self.hyperresistivity_Lambda_r = hyp['Lambda']['r']
                self.hyperresistivity_Lambda_t = hyp['Lambda']['t']
            elif 'Lambda0' in hyp:
                self.hyperresistivity_Lambda0 = float(hyp['Lambda0'])
                self.hyperresistivity_grad_j_tot_max = float(hyp['grad_j_tot_max'])
                self.hyperresistivity_min_duration = float(hyp['min_duration'])


    def setHyperresistivity(self, Lambda, radius=None, times=None):
        """
        Enable the hyperresistive diffusion term and specify the
        transport coefficient ``Lambda``.

        :param Lambda: Diffusion coefficient.
        :param radius: Radial grid on which  ``Lambda`` is specified (if any).
        :param times:  Time grid on which ``Lambda`` is specified (if any).
        """
        d, r, t = self._setPrescribedData(data=Lambda, radius=radius, times=times)

        self.hyperresistivity_type = HYPERRESISTIVITY_TYPE_PRESCRIBED
        self.hyperresistivity_Lambda_x = d
        self.hyperresistivity_Lambda_r = r
        self.hyperresistivity_Lambda_t = t


    def setHyperresistivityAdaptive(
        self, grad_j_tot_max, Lambda0, min_duration=0.5e-3
    ):
        """
        Enable the adaptive hyperresistive diffusion term, which triggers when
        the current density gradient grows above a given threshold locally. The
        term is applied until the current density drops below the threshold,
        and at least for ``min_duration`` seconds.

        :param grad_j_tot_max: Maximum current density gradient which must be exceeded for the term to be triggered.
        :param Lambda0:        Value of hyperresistivity to apply.
        :param min_duration:   Minimum duration of the hyperresistive term (in seconds).
        """
        self.hyperresistivity_type = HYPERRESISTIVITY_TYPE_ADAPTIVE
        self.hyperresistivity_grad_j_tot_max = grad_j_tot_max
        self.hyperresistivity_Lambda0 = Lambda0
        self.hyperresistivity_min_duration = min_duration


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this PoloidalFlux object.
        """
        data = {
            'hyperresistivity': {
                'type': self.hyperresistivity_type
            }
        }

        if self.hyperresistivity_type == HYPERRESISTIVITY_TYPE_PRESCRIBED:
            data['hyperresistivity']['Lambda'] = {
                'x': self.hyperresistivity_Lambda_x,
                'r': self.hyperresistivity_Lambda_r,
                't': self.hyperresistivity_Lambda_t
            }
        elif self.hyperresistivity_type == HYPERRESISTIVITY_TYPE_ADAPTIVE:
            data['hyperresistivity']['Lambda0'] = self.hyperresistivity_Lambda0
            data['hyperresistivity']['grad_j_tot_max'] = self.hyperresistivity_grad_j_tot_max
            data['hyperresistivity']['min_duration'] = self.hyperresistivity_min_duration

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.hyperresistivity_type == HYPERRESISTIVITY_TYPE_NEGLECT:
            pass
        elif self.hyperresistivity_type == HYPERRESISTIVITY_TYPE_PRESCRIBED:
            self._verifySettingsPrescribedData('psi_p hyperresistivity', data=self.hyperresistivity_Lambda_x, radius=self.hyperresistivity_Lambda_r, times=self.hyperresistivity_Lambda_t)
        elif self.hyperresistivity_type == HYPERRESISTIVITY_TYPE_ADAPTIVE:
            if not np.isscalar(self.hyperresistivity_Lambda0):
                raise EquationException(f"The hyperresistivity parameter 'Lambda0' must be a scalar. Current value: {self.hyperresistivity_Lambda0}.")
            if not np.isscalar(self.hyperresistivity_grad_j_tot_max):
                raise EquationException(f"The hyperresistivity parameter 'Lambda0' must be a scalar. Current value: {self.hyperresistivity_grad_j_tot_max}.")
            if not np.isscalar(self.hyperresistivity_min_duration):
                raise EquationException(f"The hyperresistivity parameter 'Lambda0' must be a scalar. Current value: {self.hyperresistivity_min_duration}.")
        else:
            raise EquationException(f"Invalid option for hyperresistivity type: {self.hyperresistivity_type}.")


