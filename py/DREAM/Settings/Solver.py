#
# Solver settings object
#################################

import numpy as np
from .. DREAMException import DREAMException
from . ToleranceSettings import ToleranceSettings
from . Preconditioner import Preconditioner


LINEAR_IMPLICIT = 1
NONLINEAR       = 2

LINEAR_SOLVER_LU      = 1
LINEAR_SOLVER_MUMPS   = 2
LINEAR_SOLVER_MKL     = 3
LINEAR_SOLVER_SUPERLU = 4
LINEAR_SOLVER_GMRES   = 5


class Solver:
    

    def __init__(self, ttype=LINEAR_IMPLICIT, linsolv=LINEAR_SOLVER_LU, maxiter=100, verbose=False):
        """
        Constructor.
        """
        self.setType(ttype)

        self.debug_printmatrixinfo = False
        self.debug_printjacobianinfo = False
        self.debug_savejacobian = False
        self.debug_savesolution = False
        self.debug_savematrix = False
        self.debug_savenumericaljacobian = False
        self.debug_saverhs = False
        self.debug_saveresidual = False
        self.debug_savesystem = False
        self.debug_timestep = 0
        self.debug_iteration = 1
        self.debug_rescaled = False

        self.backupsolver = None
        self.tolerance = ToleranceSettings()
        self.preconditioner = Preconditioner()
        self.setOption(linsolv=linsolv, maxiter=maxiter, verbose=verbose)


    def setDebug(self, printmatrixinfo=False, printjacobianinfo=False, savejacobian=False,
                 savesolution=False, savematrix=False, savenumericaljacobian=False, saverhs=False,
                 saveresidual=False, savesystem=False, rescaled=False, timestep=0, iteration=1):
        """
        Enable output of debug information.

        :param int timestep:   Index of time step to generate debug info for. If ``0``, debug info is generated in every (iteration of every) time step.
        :param int savesystem: Save full equation system as a DREAMOutput file in the most recent iteration/time step.

        LINEAR SOLVER
        :param bool printmatrixinfo: If ``True``, calls ``PrintInfo()`` on the linear operator matrix.
        :param bool savematrix:      If ``True``, saves the linear operator matrix using a PETSc viewer.
        :param bool saverhs:         If ``True``, saves the right-hand side vector to a ``.mat`` file.

        NON-LINEAR SOLVER
        :param bool printjacobianinfo:     If ``True``, calls ``PrintInfo()`` on the jacobian matrix.
        :param bool savejacobian:          If ``True``, saves the jacobian matrix using a PETSc viewer.
        :param bool savesolution:          If ``True``, saves the solution vector to a ``.mat`` file.
        :param bool savenumericaljacobian: If ``True``, evaluates the jacobian matrix numerically and saves it using a PETSc viewer.
        :param bool saveresidual:          If ``True``, saves the residual vector to a ``.mat`` file.
        :param bool rescaled:              If ``True``, saves the rescaled versions of the jacobian/solution/residual.
        :param int iteration:              Index of iteration to save debug info for. If ``0``, saves in all iterations. If ``timestep`` is ``0``, this parameter is always ignored.
        """
        self.debug_printmatrixinfo = printmatrixinfo
        self.debug_printjacobianinfo = printjacobianinfo
        self.debug_savejacobian = savejacobian
        self.debug_savesolution = savesolution
        self.debug_savematrix = savematrix
        self.debug_savenumericaljacobian = savenumericaljacobian
        self.debug_saverhs = saverhs
        self.debug_saveresidual = saveresidual
        self.debug_savesystem = savesystem
        self.debug_rescaled = rescaled
        self.debug_timestep = timestep
        self.debug_iteration = iteration


    def setBackupSolver(self, backup):
        """
        Set the backup linear solver to use in case the main linear
        solver fails. Set to ``None`` to disable (default).
        """
        self.backupsolver = backup


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
        """
        Specifies which type of solver to use (either ``LINEAR_IMPLICIT``
        or ``NONLINEAR``).
        """
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
        
        if 'maxiter' in data:
            self.maxiter = int(data['maxiter'])

        if 'verbose' in data:
            self.verbose = bool(data['verbose'])

        if 'tolerance' in data:
            self.tolerance.fromdict(data['tolerance'])

        if 'preconditioner' in data:
            self.preconditioner.fromdict(data['preconditioner'])

        if 'backupsolver' in data:
            self.backupsolver = int(data['backupsolver'])

        if 'debug' in data:
            flags = ['printmatrixinfo', 'printjacobianinfo', 'savejacobian', 'savesolution', 'savematrix', 'savenumericaljacobian', 'saverhs', 'saveresidual', 'savesystem', 'rescaled']

            for f in flags:
                if f in data['debug']:
                    setattr(self, 'debug_{}'.format(f), bool(data['debug'][f]))

            if 'timestep' in data['debug']:
                self.debug_timestep = int(data['debug']['timestep'])
            if 'iteration' in data['debug']:
                self.debug_iteration = int(data['debug']['iteration'])

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

        data['preconditioner'] = self.preconditioner.todict()

        if self.type == LINEAR_IMPLICIT:
            data['debug'] = {
                'printmatrixinfo': self.debug_printmatrixinfo,
                'savematrix': self.debug_savematrix,
                'saverhs': self.debug_saverhs,
                'savesystem': self.debug_savesystem,
                'timestep': self.debug_timestep
            }
        elif self.type == NONLINEAR:
            data['tolerance'] = self.tolerance.todict()
            data['debug'] = {
                'printjacobianinfo': self.debug_printjacobianinfo,
                'savejacobian': self.debug_savejacobian,
                'savesolution': self.debug_savesolution,
                'savenumericaljacobian': self.debug_savenumericaljacobian,
                'saveresidual': self.debug_saveresidual,
                'savesystem': self.debug_savesystem,
                'rescaled': self.debug_rescaled,
                'timestep': self.debug_timestep,
                'iteration': self.debug_iteration
            }

            if self.backupsolver is not None:
                data['backupsolver'] = self.backupsolver

        return data


    def verifySettings(self):
        """
        Verifies that the settings of this object are consistent.
        """
        if self.type == LINEAR_IMPLICIT:
            self.verifyLinearSolverSettings()

            if type(self.debug_printmatrixinfo) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_printmatrixinfo': {}. Expected boolean.".format(type(self.debug_printmatrixinfo)))
            elif type(self.debug_savematrix) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_savematrix': {}. Expected boolean.".format(type(self.debug_savematrix)))
            elif type(self.debug_saverhs) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_saverhs': {}. Expected boolean.".format(type(self.debug_saverhs)))
            elif type(self.debug_timestep) != int:
                raise DREAMException("Solver: Invalid type of parameter 'debug_timestep': {}. Expected integer.".format(type(self.debug_timestep)))

        elif self.type == NONLINEAR:
            if type(self.maxiter) != int:
                raise DREAMException("Solver: Invalid type of parameter 'maxiter': {}. Expected integer.".format(type(self.maxiter)))
            elif type(self.verbose) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'verbose': {}. Expected boolean.".format(type(self.verbose)))

            if type(self.debug_printjacobianinfo) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_printjacobianinfo': {}. Expected boolean.".format(type(self.debug_printjacobianinfo)))
            elif type(self.debug_savejacobian) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_savejacobian': {}. Expected boolean.".format(type(self.debug_savejacobian)))
            elif type(self.debug_savesolution) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_savesolution': {}. Expected boolean.".format(type(self.debug_savesolution)))
            elif type(self.debug_saverhs) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_saverhs': {}. Expected boolean.".format(type(self.debug_saverhs)))
            elif type(self.debug_saveresidual) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_saveresidual': {}. Expected boolean.".format(type(self.debug_saveresidual)))
            elif type(self.debug_rescaled) != bool:
                raise DREAMException("Solver: Invalid type of parameter 'debug_rescaled': {}. Expected boolean.".format(type(self.debug_rescaled)))
            elif type(self.debug_timestep) != int:
                raise DREAMException("Solver: Invalid type of parameter 'debug_timestep': {}. Expected integer.".format(type(self.debug_timestep)))
            elif type(self.debug_iteration) != int:
                raise DREAMException("Solver: Invalid type of parameter 'debug_iteration': {}. Expected boolean.".format(type(self.debug_iteration)))

            self.tolerance.verifySettings()
            self.verifyLinearSolverSettings()
        else:
            raise DREAMException("Solver: Unrecognized solver type: {}.".format(self.type))

        self.preconditioner.verifySettings()


    def verifyLinearSolverSettings(self):
        """
        Verifies the settings for the linear solver (which is used
        by both the 'LINEAR_IMPLICIT' and 'NONLINEAR' solvers).
        """
        solv = [LINEAR_SOLVER_LU, LINEAR_SOLVER_MUMPS, LINEAR_SOLVER_MKL, LINEAR_SOLVER_SUPERLU, LINEAR_SOLVER_GMRES]
        if self.linsolv not in solv:
            raise DREAMException("Solver: Unrecognized linear solver type: {}.".format(self.linsolv))
        elif self.backupsolver is not None and self.backupsolver not in solv:
            raise DREAMException("Solver: Unrecognized backup linear solver type: {}.".format(self.backupsolver))


