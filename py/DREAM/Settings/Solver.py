#
# Solver settings object
#################################

import numpy as np
from .. DREAMException import DREAMException
from . ToleranceSettings import ToleranceSettings


LINEAR_IMPLICIT = 1
NONLINEAR       = 2

LINEAR_SOLVER_LU    = 1
LINEAR_SOLVER_MUMPS = 2


class Solver:
    

    def __init__(self, ttype=LINEAR_IMPLICIT, linsolv=LINEAR_SOLVER_LU, maxiter=100, verbose=False):
        """
        Constructor.
        """
        self.setType(ttype)

        self.tolerance = ToleranceSettings()
        self.setOption(linsolv=linsolv, maxiter=maxiter, verbose=verbose)


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
        print("WARNING: The 'Solver.setTolerance()' method is deprecated. Please use 'Solver.tolerance.set(reltol=...)' instead.")
        self.tolerance.set(reltol=reltol)


    def setVerbose(self, verbose):
        """
        If 'True', generates excessive output during nonlinear solve.
        """
        self.setOption(verbose=verbose)


    def setOption(self, linsolv=None, maxiter=None, verbose=None):
        """
        Sets a solver option.
        """
        if linsolv is not None:
            self.linsolv = linsolv
        if maxiter is not None:
            self.maxiter = maxiter
        if verbose is not None:
            self.verbose = verbose

        self.verifySettings()


    def setType(self, ttype):
        if ttype == LINEAR_IMPLICIT:
            self.type = ttype
        elif ttype == NONLINEAR:
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
        self.verbose = bool(data['verbose'])

        if 'tolerance' in data:
            self.tolerance.fromdict(data['tolerance'])

        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this Solver object.
        """
        if verify:
            self.verifySettings()

        data = {
            'type': self.type,
            'linsolv': self.linsolv,
            'maxiter': self.maxiter,
            'verbose': self.verbose
        }

        if self.type == NONLINEAR:
            data['tolerance'] = self.tolerance.todict()

        return data


    def verifySettings(self):
        """
        Verifies that the settings of this object are consistent.
        """
        if self.type == LINEAR_IMPLICIT:
            self.verifyLinearSolverSettings()
        elif self.type == NONLINEAR:
            if type(self.maxiter) != int:
                raise DREAMException("Solver: Invalid type of parameter 'maxiter': {}. Expected integer.".format(self.maxiter))
            elif type(self.verbose) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'verbose': {}. Expected boolean.".format(self.verbose))

            self.tolerance.verifySettings()
            self.verifyLinearSolverSettings()
        else:
            raise DREAMException("Solver: Unrecognized solver type: {}.".format(self.type))


    def verifyLinearSolverSettings(self):
        """
        Verifies the settings for the linear solver (which is used
        by both the 'LINEAR_IMPLICIT' and 'NONLINEAR' solvers).
        """
        solv = [LINEAR_SOLVER_LU, LINEAR_SOLVER_MUMPS]
        if self.linsolv not in solv:
            raise DREAMException("Solver: Unrecognized linear solver type: {}.".format(self.linsolv))


