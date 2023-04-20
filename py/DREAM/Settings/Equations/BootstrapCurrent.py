# Settings for the Bootstrap current
#S
#
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity

BOOTSTRAP_MODE_DISABLED = 1
BOOTSTRAP_MODE_ENABLED = 2

BOOTSTRAP_BC_BACKWARDS = 1
BOOTSTRAP_BC_ZERO = 2


class BootstrapCurrent(UnknownQuantity):

    def __init__(self, settings, mode=BOOTSTRAP_MODE_DISABLED, bc=BOOTSTRAP_BC_BACKWARDS):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.setMode(mode)
        self.setBoundaryCondition(bc)


    def setMode(self, mode):
        """
        Specifies whether to include the bootstrap current as a contribution to the total
        current density. If enabled, this contribution is calculated using the Redl-Sauter
        model, which is based on A. Redl et al (DOI: https://doi.org/10.1063/5.0012664).
        """
        if mode in [BOOTSTRAP_MODE_DISABLED, BOOTSTRAP_MODE_ENABLED]:
            self.mode = mode
        else:
            raise EquationException("j_bs: Unrecognized bootstrap current mode: {}".format(mode))

    def setBoundaryCondition(self, bc):
        """
        Specifies whether to set the boundary condition using a forward Euler method at the edge,
        or to assume that densities/temperatures vanish outside the edge and use central difference.
        """
        if bc in [BOOTSTRAP_BC_BACKWARDS, BOOTSTRAP_BC_ZERO]:
            self.bc = bc
        else:
            raise EquationException("j_bs: Unrecognized bootstrap current boundary condition: {}".format(mode))


    def fromdict(self, data):
        """
        Set options from a dictionary.
        """
        if 'mode' in data:
            self.setMode(data['mode'])
        if 'bc' in data:
            self.setBoundaryCondition(data['bc'])

    def todict(self):
        """
        Returns a Python dictionary containing the settings of this BootstrapCurrent object.
        """
        data = {
            'mode': self.mode,
            'bc':   self.bc
        }
        return data

    def verifySettings(self):
        pass
