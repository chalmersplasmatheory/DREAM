#
# Solver settings object
#################################

import numpy as np


LINEAR_IMPLICIT = 1
NONLINEAR_SNES  = 2


class Solver:
    

    def __init__(self, ttype=LINEAR_IMPLICIT, maxiter=100, reltol=1e-6, verbose=False):
        """
        Constructor.
        """
        self.setType(ttype)

        self.setOption(maxiter=maxiter, reltol=reltol, verbose=verbose)


    def setMaxIterations(self, maxiter):
        """
        Set maximum number of allowed nonlinear iterations.
        """
        self.setOption(maxiter=maxiter)


    def setTolerance(self, reltol):
        """
        Set relative tolerance for nonlinear solve.
        """
        self.setOption(reltol=reltol)


    def setVerbose(self, verbose):
        """
        If 'True', generates excessive output during nonlinear solve.
        """
        self.setOption(verbose=verbose)


    def setOption(self, maxiter=None, reltol=None, verbose=None):
        """
        Sets a solver option.
        """
        if maxiter is not None:
            self.maxiter = maxiter
        if reltol is not None:
            self.reltol = reltol
        if verbose is not None:
            self.verbose = verbose

        self.verifySettings()


    def setType(self, ttype):
        if ttype == LINEAR_IMPLICIT:
            self.type = ttype
        elif ttype == NONLINEAR_SNES:
            self.type = ttype
        else:
            raise DREAMException("Solver: Unrecognized solver type: {}.".format(ttype))


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
            'type': self.type,
            'maxiter': self.maxiter,
            'reltol': self.reltol,
            'verbose': self.verbose
        }


    def verifySettings(self):
        """
        Verifies that the settings of this object are consistent.
        """
        if self.type == LINEAR_IMPLICIT: pass
        elif self.type == NONLINEAR_SNES:
            if type(self.maxiter) != int:
                raise DREAMException("Solver: Invalid type of parameter 'maxiter': {}. Expected integer.".format(self.maxiter))
            elif type(self.reltol) != float and type(self.reltol) != int:
                raise DREAMException("Solver: Invalid type of parameter 'reltol': {}. Expected float.".format(self.reltol))
            elif type(self.verbose) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'verbose': {}. Expected boolean.".format(self.verbose))
        else:
            raise DREAMException("Solver: Unrecognized solver type: {}.".format(ttype))


