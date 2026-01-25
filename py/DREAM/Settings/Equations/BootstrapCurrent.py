# Settings for the Bootstrap current
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity

from . PrescribedInitialParameter import PrescribedInitialParameter

BOOTSTRAP_MODE_DISABLED = 1
BOOTSTRAP_MODE_REDL = 2

BOOTSTRAP_INIT_MODE_TOTAL = 1 
BOOTSTRAP_INIT_MODE_OHMIC = 2

class BootstrapCurrent(UnknownQuantity):

    def __init__(self, settings, mode=BOOTSTRAP_MODE_DISABLED, initMode=BOOTSTRAP_INIT_MODE_TOTAL):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.mode = None
        self.initMode = None

        self.setMode(mode)
        self.setInitMode(initMode)


    def setMode(self, mode, initMode=None):
        """
        Specifies whether to include the bootstrap current as a contribution to the total
        current density. If enabled, this contribution is calculated using the Redl-Sauter
        model, which is based on A. Redl et al (DOI: https://doi.org/10.1063/5.0012664).
        """
        if mode is True:
            self.mode = BOOTSTRAP_REDL
        elif mode is False:
            self.mode = BOOTSTRAP_DISABLED
        elif mode in [BOOTSTRAP_MODE_DISABLED, BOOTSTRAP_MODE_REDL]:
            self.mode = int(mode)
        else:
            print(type(mode), mode)
            raise EquationException(f"j_bs: Unrecognized bootstrap current mode: {mode}")
        
        if initMode is not None:
            self.setInitMode(initMode)


    def setInitMode(self, initMode):
        """
        Specifies whether the initial/prescribed total current density is given as only the
        Ohmic current (so the bootstrap current is added to obtain the correct total current)
        or if it is given as the sum of the Ohmic and bootstrap current (so the bootstrap
        current is subtracted from the total current when calculating the corresponding electric
        field).
        """
        if initMode in [BOOTSTRAP_INIT_MODE_OHMIC, BOOTSTRAP_INIT_MODE_TOTAL]:
            self.initMode = int(initMode)
        else:
            raise EquationException(f"j_bs: Unrecognized bootstrap current initialization mode: {initMode}")


    def fromdict(self, data):
        """
        Set options from a dictionary.
        """
        if 'mode' in data:
            self.setMode(data['mode'])
        if 'init_mode' in data:
            self.setInitMode(data['init_mode'])


    def todict(self):
        """
        Returns a Python dictionary containing the settings of this BootstrapCurrent object.
        """
        data = {'mode': self.mode, 'init_mode': self.initMode}
        return data


    def verifySettings(self):
        """
        Verify that the settings concerning the bootstrap current are correctly set.
        """
        if self.mode == BOOTSTRAP_MODE_DISABLED and self.initMode == BOOTSTRAP_INIT_MODE_OHMIC:
            print("WARNING: Bootstrap current is disabled, but its initialization mode was adjusted!")
