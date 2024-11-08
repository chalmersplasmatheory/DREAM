import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings


HYPERRESISTIVITY_MODE_NEGLECT = 1
HYPERRESISTIVITY_MODE_PRESCRIBED = 2
HYPERRESISTIVITY_MODE_ADAPTIVE = 3


class PoloidalFlux(UnknownQuantity,PrescribedParameter):
    

    def __init__(self, settings):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.hyperresistivity_mode = HYPERRESISTIVITY_MODE_NEGLECT
        self.hyperresistivity_Lambda_x = None
        self.hyperresistivity_Lambda_r = None
        self.hyperresistivity_Lambda_t = None

        self.hyperresistivity_grad_j_tot_max = None
        self.hyperresistivity_gradient_normalized = False
        self.hyperresistivity_Lambda0 = None
        self.hyperresistivity_min_duration = None


    def setHyperresistivity(self, Lambda, radius=None, times=None):
        """
        Enable the hyperresistive diffusion term and specify the
        transport coefficient ``Lambda``.

        :param Lambda: Diffusion coefficient.
        :param radius: Radial grid on which  ``Lambda`` is specified (if any).
        :param times:  Time grid on which ``Lambda`` is specified (if any).
        """
        d, r, t = self._setPrescribedData(data=Lambda, radius=radius, times=times)

        self.hyperresistivity_mode = HYPERRESISTIVITY_MODE_PRESCRIBED
        self.hyperresistivity_Lambda_x = d
        self.hyperresistivity_Lambda_r = r
        self.hyperresistivity_Lambda_t = t


    def setHyperresistivityAdaptive(
        self, Lambda0, grad_j_tot_max=None,
        grad_j_tot_max_norm=None, min_duration=0.5e-3
    ):
        """
        Enable the adaptive hyperresistive diffusion term, which triggers when
        the current density gradient grows above a given threshold locally. The
        term is applied until the current density drops below the threshold,
        and at least for ``min_duration`` seconds.

        :param Lambda0:             Value of hyperresistivity to apply.
        :param grad_j_tot_max:      Maximum current density gradient which must be exceeded for the term to be triggered.
        :param grad_j_tot_max_norm: Maximum current density gradient (normalized to average current density) which must be exceeded.
        :param min_duration:        Minimum duration of the hyperresistive term (in seconds).
        """
        self.hyperresistivity_mode = HYPERRESISTIVITY_MODE_ADAPTIVE

        if grad_j_tot_max:
            self.hyperresistivity_grad_j_tot_max = grad_j_tot_max
            self.hyperresistivity_gradient_normalized = False
        elif grad_j_tot_max_norm:
            self.hyperresistivity_grad_j_tot_max = grad_j_tot_max_norm
            self.hyperresistivity_gradient_normalized = True
        else:
            raise EquationException("One of 'grad_j_tot_max' and 'grad_j_tot_max_norm' must be specified.")

        self.hyperresistivity_Lambda0 = Lambda0
        self.hyperresistivity_min_duration = min_duration


    def fromdict(self, data):
        """
        Set all options from a dictionary.
        """
        if 'hyperresistivity' in data:
            hyp = data['hyperresistivity']
            self.hyperresistivity_mode = int(hyp['mode'])

            if 'Lambda' in hyp:
                self.hyperresistivity_Lambda_x = hyp['Lambda']['x']
                self.hyperresistivity_Lambda_r = hyp['Lambda']['r']
                self.hyperresistivity_Lambda_t = hyp['Lambda']['t']
            elif 'Lambda0' in hyp:
                self.hyperresistivity_Lambda0 = float(hyp['Lambda0'])
                self.hyperresistivity_grad_j_tot_max = float(hyp['grad_j_tot_max'])
                self.hyperresistivity_gradient_normalized = bool(hyp['gradient_normalized'])
                self.hyperresistivity_min_duration = float(hyp['min_duration'])


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this PoloidalFlux object.
        """
        hypres = { 'mode': self.hyperresistivity_mode }

        if self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_PRESCRIBED:
            hypres['Lambda'] = {
                'x': self.hyperresistivity_Lambda_x,
                'r': self.hyperresistivity_Lambda_r,
                't': self.hyperresistivity_Lambda_t
            }
        elif self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_ADAPTIVE:
            hypres['Lambda0'] = self.hyperresistivity_Lambda0
            hypres['grad_j_tot_max'] = self.hyperresistivity_grad_j_tot_max
            hypres['gradient_normalized'] = self.hyperresistivity_gradient_normalized
            hypres['min_duration'] = self.hyperresistivity_min_duration

        return { 'hyperresistivity': hypres }


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_NEGLECT:
            pass
        elif self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_PRESCRIBED:
            self._verifySettingsPrescribedData('psi_p hyperresistivity', data=self.hyperresistivity_Lambda_x, radius=self.hyperresistivity_Lambda_r, times=self.hyperresistivity_Lambda_t)
        elif self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_ADAPTIVE:
            if not np.isscalar(self.hyperresistivity_Lambda0):
                raise EquationException(f"The hyperresistivity parameter 'Lambda0' must be a scalar. Current value: {self.hyperresistivity_Lambda0}.")
            if not np.isscalar(self.hyperresistivity_grad_j_tot_max):
                raise EquationException(f"The hyperresistivity parameter 'grad_j_tot_max' must be a scalar. Current value: {self.hyperresistivity_grad_j_tot_max}.")
            if not np.isscalar(self.hyperresistivity_min_duration):
                raise EquationException(f"The hyperresistivity parameter 'min_duration' must be a scalar. Current value: {self.hyperresistivity_min_duration}.")
        else:
            raise EquationException(f"Invalid option for hyperresistivity mode: {self.hyperresistivity_mode}.")


