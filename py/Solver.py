#
# Solver settings object
#################################

import numpy as np


TYPE_LINEAR_IMPLICIT = 1
TYPE_NONLINEAR_SNES  = 2


class Solver:
    

    TYPE_LINEAR_IMPLICIT = 1
    TYPE_NONLINEAR_SNES  = 2


    def __init__(self, ttype=1):
        """
        Constructor.
        """
        self.setType(ttype)


    def setType(self, ttype):
        if ttype == self.TYPE_LINEAR_IMPLICIT:
            self.type = ttype
        elif ttype == self.TYPE_NONLINEAR_SNES:
            self.type = ttype
        else:
            raise DREAMException("Solve: Unrecognized solver type: {}.".format(ttype))


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
        if self.type == self.TYPE_LINEAR_IMPLICIT: pass
        elif self.type == self.TYPE_NONLINEAR_SNES: pass
        else:
            raise DREAMException("Solve: Unrecognized solver type: {}.".format(ttype))


