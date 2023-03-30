# Settings for the Bootstrap current
#
# NOTE: move this to e.g. ElectricField/OhmicCurrent??
#
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity

BOOTSTRAP_MODE_DISABLED = 1
BOOTSTRAP_MODE_ENABLED = 2

class BootstrapCurrent(UnknownQuantity):

    def __init__(self, settings, mode=BOOTSTRAP_MODE_DISABLED):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.setMode(mode)


    def setMode(self, mode):
        """
        Specifies whether to include the bootstrap current as a contribution to the total
        current density. If enabled, this contribution is calculated using the Redl-Sauter
        model, which is based on A. Redl et al (DOI: https://doi.org/10.1063/5.0012664).
        """
        if mode in [BOOTSTRAP_MODE_DISABLED, BOOTSTRAP_MODE_DISABLED]:
            self.mode = mode
        else:
            raise EquationException("j_bs: Unrecognized bootstrap current mode: {}".format(mode))


    def fromdict(self, data):
        """
        Set options from a dictionary.
        """
        if 'mode' in data:
            self.setMode(mode)

    def todict(self):
        """
        Returns a Python dictionary containing the settings of this BootstrapCurrent object.
        """
        data = { 'mode': self.mode }
        return data

    def verifySettings(self):
        pass
