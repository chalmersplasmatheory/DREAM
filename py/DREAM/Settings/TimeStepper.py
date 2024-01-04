#
# TimeStepper settings object
# #################################

import numpy as np

from .. DREAMException import DREAMException
from . ToleranceSettings import ToleranceSettings


TYPE_CONSTANT = 1
TYPE_ADAPTIVE = 2
TYPE_IONIZATION = 3


class TimeStepper:
    
    def __init__(self, ttype=1, checkevery=0, tmax=None, dt=None, nt=None, nSaveSteps=0, reltol=1e-2, verbose=False, constantstep=False, terminatefunc=None):
        """
        Constructor.
        """
        self.set(ttype=ttype, checkevery=checkevery, tmax=tmax, dt=dt, nt=nt, nSaveSteps=nSaveSteps, reltol=reltol, verbose=verbose, constantstep=constantstep, terminatefunc=None)
        

    def set(self, ttype=1, checkevery=0, tmax=None, dt=None, nt=None, nSaveSteps=0, reltol=1e-2, verbose=False, constantstep=False, minsavedt=0, terminatefunc=None):
        """
        Set properties of the time stepper.
        """
        self.type = int(ttype)

        self.setCheckInterval(checkevery)
        self.setTmax(tmax)
        self.setDt(dt)
        self.setNt(nt)
        self.setMinSaveTimestep(minsavedt)
        self.setNumberOfSaveSteps(nSaveSteps)
        self.setVerbose(verbose)
        self.setConstantStep(constantstep)       
        self.tolerance = ToleranceSettings()
        self.tolerance.set(reltol=reltol)
        self.terminatefunc = terminatefunc
        
        self.dtmax = None
        self.automaticstep = None
        self.safetyfactor = None


    def __contains__(self, item):
        return (item in self.todict(False))


    def __getitem__(self, key):
        return self.todict(False)[key]


    ######################
    # SETTERS
    ######################
    def setCheckInterval(self, checkevery):
        if checkevery < 0:
            raise DREAMException("TimeStepper: Invalid value assigned to 'checkevery': {}".format(checkevery))
        
        self.checkevery = int(checkevery)


    def setConstantStep(self, constantstep):
        self.constantstep = bool(constantstep)


    def setDt(self, dt):
        if dt is None:
            self.dt = None
            return

        if dt < 0 or (dt == 0 and self.type != TYPE_IONIZATION):
            raise DREAMException("TimeStepper: Invalid value assigned to 'dt': {}".format(tmax))
        if self.nt is not None and dt > 0:
            raise DREAMException("TimeStepper: 'dt' may not be set alongside 'nt'.")
            
        self.dt = float(dt)


    def setMinSaveTimestep(self, dt):
        """
        For the adapative ionization-based time stepper, sets the minimum
        time which must elapse between two saved time steps.
        """
        self.minsavedt = dt


    def setNt(self, nt):
        if nt is None:
            self.nt = None
            return

        if nt <= 0:
            raise DREAMException("TimeStepper: Invalid value assigned to 'dt': {}".format(tmax))
        if self.dt is not None and self.dt > 0:
            raise DREAMException("TimeStepper: 'nt' may not be set alongside 'dt'.")
            
        self.nt = int(nt)


    def setNumberOfSaveSteps(self, nSaveSteps):
        """
        Sets the number of time steps to save to the output file.
        This number must be <= Nt. If 0, all time steps are saved.
        """
        self.nSaveSteps = nSaveSteps


    def setRelTol(self, reltol): self.setRelativeTolerance(reltol=reltol)


    def setRelativeTolerance(self, reltol):
        if reltol <= 0:
            raise DREAMException("TimeStepper: Invalid value assigned to 'reltol': {}".format(reltol))

        self.tolerance.set(reltol=float(reltol))


    def setTerminationFunction(self, func):
        """
        Sets the Python function to call in order to determine when terminate
        the time stepping. **NOTE**: This functionality is only available when
        DREAM is compiled and run as a Python library.

        :param func: Python function determining when to terminate time stepping. Takes a libdreampyface 'Simulation' object as input and returns a bool.
        """
        self.terminatefunc = func


    def setTmax(self, tmax):
        if tmax is None:
            self.tmax = None
            return

        if tmax <= 0:
            raise DREAMException("TimeStepper: Invalid value assigned to 'tmax': {}".format(tmax))

        self.tmax = float(tmax)


    def setType(self, ttype, *args, **kwargs):
        if ttype not in [TYPE_CONSTANT, TYPE_ADAPTIVE, TYPE_IONIZATION]:
            raise DREAMException("TimeStepper: Unrecognized time stepper type specified: {}".format(ttype))

        if ttype in [TYPE_ADAPTIVE, TYPE_IONIZATION]:
            self.nt = None

        self.type = int(ttype)
        
        if ttype == TYPE_IONIZATION:
            self.setIonization(*args, **kwargs)
    

    def setIonization(self, dt0=0, dtmax=0, tmax=None, automaticstep=1e-12, safetyfactor=50):
        """
        Select and set parameters for the ionization time stepper.
        """
        self.type = TYPE_IONIZATION
        self.dt = dt0
        self.dtmax = dtmax
        self.automaticstep = automaticstep
        self.safetyfactor = safetyfactor

        if tmax is not None:
            self.tmax = tmax


    def setVerbose(self, verbose=True):
        self.verbose = bool(verbose)


    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        def scal(v):
            if type(v) == np.ndarray: return v[0]
            else: return v

        self.type = data['type']
        self.tmax = data['tmax']

        if type(self.type) == np.ndarray: self.type = int(self.type.flatten()[0])
        if type(self.tmax) == np.ndarray: self.tmax = float(self.tmax.flatten()[0])

        if 'automaticstep' in data: self.automaticstep = float(scal(data['automaticstep']))
        if 'checkevery' in data: self.checkevery = int(scal(data['checkevery']))
        if 'constantstep' in data: self.constantstep = bool(scal(data['constantstep']))
        if 'dt' in data: self.dt = float(scal(data['dt']))
        if 'dtmax' in data: self.dtmax = float(scal(data['dtmax']))
        if 'minsavedt' in data: self.minsavedt = float(scal(data['minsavedt']))
        if 'nt' in data: self.nt = int(scal(data['nt']))
        if 'nsavesteps' in data: self.nSaveSteps = int(scal(data['nsavesteps']))
        if 'verbose' in data: self.verbose = bool(scal(data['verbose']))
        if 'safetyfactor' in data: self.safetyfactor = float(scal(data['safetyfactor']))
        if 'tolerance' in data: self.tolerance.fromdict(data['tolerance'])
        if 'terminatefunc' in data: self.terminatefunc = data['terminatefunc']
        
        self.verifySettings()


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this TimeStepper object.
        """
        if verify:
            self.verifySettings()

        data = {
            'type': self.type,
            'tmax': self.tmax
        }

        if self.dt is not None: data['dt'] = self.dt

        if self.type == TYPE_CONSTANT:
            if self.nt is not None: data['nt'] = self.nt
            data['nsavesteps'] = int(self.nSaveSteps)

            if self.terminatefunc != None:
                data['terminatefunc'] = self.terminatefunc
        elif self.type == TYPE_ADAPTIVE:
            data['checkevery'] = self.checkevery
            data['constantstep'] = self.constantstep
            data['tolerance'] = self.tolerance.todict()
            data['verbose'] = self.verbose
        elif self.type == TYPE_IONIZATION:
            if self.dtmax is not None: data['dtmax'] = self.dtmax
            data['automaticstep'] = self.automaticstep
            data['safetyfactor'] = self.safetyfactor
            data['minsavedt'] = self.minsavedt

        return data


    def verifySettings(self):
        """
        Verify that the TimeStepper settings are consistent.
        """
        if self.type == TYPE_CONSTANT:
            if self.tmax is None or self.tmax <= 0:
                raise DREAMException("TimeStepper constant: 'tmax' must be set to a value > 0.")
            
            # Verify that _exactly_ one of 'dt' and 'nt' is
            # set to a valid value
            dtSet = (self.dt is not None and self.dt > 0)
            ntSet = (self.nt is not None and self.nt > 0)

            if dtSet and ntSet:
                raise DREAMException("TimeStepper constant: Exactly one of 'dt' and 'nt' must be > 0.")

            if self.nSaveSteps < 0 or (ntSet and self.nSaveSteps > self.nt):
                raise DREAMException("TimeStepper constant: Invalid value assigned to 'nSaveSteps'. Must between 0 and nt.")
        elif self.type == TYPE_ADAPTIVE:
            if self.tmax is None or self.tmax <= 0:
                raise DREAMException("TimeStepper adaptive: 'tmax' must be set to a value > 0.")
            elif self.nt is not None:
                raise DREAMException("TimeStepper adaptive: 'nt' cannot be used with the adaptive time stepper.")

            if type(self.checkevery) != int or self.checkevery < 0:
                raise DREAMException("TimeStepper adaptive: 'checkevery' must be a non-negative integer.")
            elif type(self.verbose) != bool:
                raise DREAMException("TimeStepper adaptive: 'verbose' must be a boolean.")
            elif type(self.constantstep) != bool:
                raise DREAMException("TimeStepper adaptive: 'constantstep' must be a boolean.")
            self.tolerance.verifySettings()
        elif self.type == TYPE_IONIZATION:
            if self.tmax is None or self.tmax <= 0:
                raise DREAMException("TimeStepper ionization: 'tmax' must be set to a value > 0.")
            elif self.dt is None or self.dt < 0:
                raise DREAMException("TimeStepper ionization: 'dt' must be set to a non-negative value.")
            elif self.dtmax is None or self.dtmax < 0:
                raise DREAMException("TimeStepper ionization: 'dtmax' must be set to a non-negative value.")
            elif self.minsavedt < 0:
                raise DREAMException("TimeStepper ionization: 'minsavedt' must be non-negative.")
        else:
            raise DREAMException("Unrecognized time stepper type selected: {}.".format(self.type))


