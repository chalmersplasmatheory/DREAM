#
# Solver settings object
#################################

import numpy as np


LINEAR_IMPLICIT = 1
NONLINEAR_SNES  = 2


class Solver:
    

    def __init__(self, ttype=LINEAR_IMPLICIT):
        """
        Constructor.
        """
        self.setType(ttype)


    def setType(self, ttype):
        if ttype == LINEAR_IMPLICIT:
            self.type = ttype
        elif ttype == NONLINEAR_SNES:
            self.type = ttype
        else:
            raise DREAMException("Solve: Unrecognized solver type: {}.".format(ttype))


    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.type = data['type']

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this Solver object.
        """
        if verify:
            self.verifySettings()

        return {
            'type': self.type
        }


    def verifySettings(self):
        """
        Verifies that the settings of this object are consistent.
        """
        if self.type == LINEAR_IMPLICIT: pass
        elif self.type == NONLINEAR_SNES: pass
        else:
            raise DREAMException("Solve: Unrecognized solver type: {}.".format(ttype))


