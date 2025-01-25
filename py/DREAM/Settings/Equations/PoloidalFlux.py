import numpy as np
from . EquationException import EquationException
from . PrescribedParameter import PrescribedParameter
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings


HYPERRESISTIVITY_MODE_NEGLECT = 1
HYPERRESISTIVITY_MODE_PRESCRIBED = 2
HYPERRESISTIVITY_MODE_ADAPTIVE = 3
HYPERRESISTIVITY_MODE_ADAPTIVE_LOCAL = 4


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
        self.hyperresistivity_dBB0 = None
        self.hyperresistivity_suppression_level = 0.9


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
        self, dBB0, grad_j_tot_max=None,
        grad_j_tot_max_norm=None,
        localized=False, suppression_level=0.9
    ):
        """
        Enable the adaptive hyperresistive diffusion term, which triggers when
        the current density gradient grows above a given threshold locally. The
        term is applied until the current density drops below the threshold,
        and at least untile the gradient reaches the 10% of the threshold.

        :param dBB0:                Value of hyperresistivity to apply.
        :param grad_j_tot_max:      Maximum current density gradient which must be exceeded for the term to be triggered.
        :param grad_j_tot_max_norm: Maximum current density gradient (normalized to average current density) which must be exceeded.
        :param localized:           Apply localized transport.
        """
        if localized:
            self.hyperresistivity_mode = HYPERRESISTIVITY_MODE_ADAPTIVE_LOCAL
        else:
            self.hyperresistivity_mode = HYPERRESISTIVITY_MODE_ADAPTIVE

        if grad_j_tot_max:
            self.hyperresistivity_grad_j_tot_max = grad_j_tot_max
            self.hyperresistivity_gradient_normalized = False
        elif grad_j_tot_max_norm:
            self.hyperresistivity_grad_j_tot_max = grad_j_tot_max_norm
            self.hyperresistivity_gradient_normalized = True
        else:
            raise EquationException("One of 'grad_j_tot_max' and 'grad_j_tot_max_norm' must be specified.")

        self.hyperresistivity_dBB0 = dBB0
        self.hyperresistivity_suppression_level = suppression_level


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
            elif 'dBB0' in hyp:
                self.hyperresistivity_dBB0 = float(hyp['dBB0'])
                self.hyperresistivity_grad_j_tot_max = float(hyp['grad_j_tot_max'])
                self.hyperresistivity_gradient_normalized = bool(hyp['gradient_normalized'])
                self.hyperresistivity_suppression_level = float(hyp['suppression_level'])


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
        elif self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_ADAPTIVE or self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_ADAPTIVE_LOCAL:
            hypres['dBB0'] = self.hyperresistivity_dBB0
            hypres['grad_j_tot_max'] = self.hyperresistivity_grad_j_tot_max
            hypres['gradient_normalized'] = self.hyperresistivity_gradient_normalized
            hypres['suppression_level'] = self.hyperresistivity_suppression_level 

        return { 'hyperresistivity': hypres }


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_NEGLECT:
            pass
        elif self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_PRESCRIBED:
            self._verifySettingsPrescribedData('psi_p hyperresistivity', data=self.hyperresistivity_Lambda_x, radius=self.hyperresistivity_Lambda_r, times=self.hyperresistivity_Lambda_t)
        elif self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_ADAPTIVE or self.hyperresistivity_mode == HYPERRESISTIVITY_MODE_ADAPTIVE_LOCAL:
            if not np.isscalar(self.hyperresistivity_dBB0):
                raise EquationException(f"The hyperresistivity parameter 'dBB0' must be a scalar. Current value: {self.hyperresistivity_dBB0}.")
            if not np.isscalar(self.hyperresistivity_grad_j_tot_max):
                raise EquationException(f"The hyperresistivity parameter 'grad_j_tot_max' must be a scalar. Current value: {self.hyperresistivity_grad_j_tot_max}.")
            if not np.isscalar(self.hyperresistivity_suppression_level):
                raise EquationException(f"The hyperresistivity parameter 'suppression_level' must be a scalar. Current value: {self.hyperresistivity_suppression_level}.")
        else:
            raise EquationException(f"Invalid option for hyperresistivity mode: {self.hyperresistivity_mode}.")


