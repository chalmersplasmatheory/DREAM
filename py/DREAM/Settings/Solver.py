#
# Solver settings object
#################################

import numpy as np
from .. DREAMException import DREAMException


LINEAR_IMPLICIT = 1
NONLINEAR       = 2
NONLINEAR_SNES  = 3

LINEAR_SOLVER_LU    = 1
LINEAR_SOLVER_GMRES = 2
LINEAR_SOLVER_MUMPS = 3


class Solver:
    

    def __init__(self, ttype=LINEAR_IMPLICIT, linsolv=LINEAR_SOLVER_LU, maxiter=100, reltol=1e-6, verbose=False):
        """
        Constructor.
        """
        self.setType(ttype)

        self.setOption(linsolv=linsolv, maxiter=maxiter, reltol=reltol, verbose=verbose)


    def setLinearSolver(self, linsolv):
        """
        Set the linear solver to use.
        """
        self.linsolv = linsolv


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

    def setOption(self, linsolv=None, maxiter=None, reltol=None, verbose=None):
        """
        Sets a solver option.
        """
        if linsolv is not None:
            self.linsolv = linsolv
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
        elif ttype == NONLINEAR:
            self.type = ttype
        elif ttype == NONLINEAR_SNES:
            self.type = ttype
        else:
            raise DREAMException("Solver: Unrecognized solver type: {}.".format(ttype))


    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        def scal(v):
            if type(v) == np.ndarray: return v[0]
            else: return v

        self.type = int(scal(data['type']))
        self.linsolv = int(data['linsolv'])
        self.maxiter = int(data['maxiter'])
        self.reltol = float(data['reltol'])
        self.verbose = bool(data['verbose'])

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
            'linsolv': self.linsolv,
            'maxiter': self.maxiter,
            'reltol': self.reltol,
            'verbose': self.verbose
        }


    def verifySettings(self):
        """
        Verifies that the settings of this object are consistent.
        """
        if self.type == LINEAR_IMPLICIT:
            self.verifyLinearSolverSettings()
        elif (self.type == NONLINEAR) or (self.type == NONLINEAR_SNES):
            if type(self.maxiter) != int:
                raise DREAMException("Solver: Invalid type of parameter 'maxiter': {}. Expected integer.".format(self.maxiter))
            elif type(self.reltol) != float and type(self.reltol) != int:
                raise DREAMException("Solver: Invalid type of parameter 'reltol': {}. Expected float.".format(self.reltol))
            elif type(self.verbose) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'verbose': {}. Expected boolean.".format(self.verbose))

            self.verifyLinearSolverSettings()
        else:
            raise DREAMException("Solver: Unrecognized solver type: {}.".format(self.type))


    def verifyLinearSolverSettings(self):
        """
        Verifies the settings for the linear solver (which is used
        by both the 'LINEAR_IMPLICIT' and 'NONLINEAR' solvers).
        """
        solv = [LINEAR_SOLVER_LU, LINEAR_SOLVER_GMRES, LINEAR_SOLVER_MUMPS]
        if self.linsolv not in solv:
            raise DREAMException("Solver: Unrecognized linear solver type: {}.".format(self.linsolv))


