# A high-level class for easily conducting
# convergence scans.
##################################################

import numpy as np
from . import DREAMIO
from . import runiface

from DREAM.DREAMSettings import DREAMSettings
from DREAM.DREAMOutput import DREAMOutput


class ConvergenceScan:


    def __init__(self, settings, inparams=None, outparams=None, scanUntilConvergence=False,
                 verbose=True):
        """
        Creates a new ConvergenceScan object with 'settings' representing
        the settings for the baseline scenario.

        :param DREAMSettings settings:    Baseline ``DREAMSettings`` object to base all convergence runs on.
        :param list inparams:             Either a string (or a list of strings), specifying the set(s) of parameters to scan, or ``None``, which sets no parameters (and they must then be explicitly set using ``addScanParameter()`` later).
        :param list outparams:            Either a string (or a list of strings), specifying the set(s) of parameters to measure for convergence. Alternatively, ``None`` clears all output parameters (which must then be explicitly set using ``addOutputParameter()`` later).
        :param bool scanUntilConvergence: If ``True``, does not limit the number of runs to do and updates the value of each parameter until the output parameter changes less than the given tolerance.
        :param bool verbose:              If ``True``, prints progress message to stdout when running.
        """
        self.settings = settings
        self.scanParameters = dict()
        self.outputParameters = dict()

        self.result = {}
        self.baselineOutput = None
        self.verbose = verbose

        # Maximum number of iterations if running scan adaptively
        self.NMAX = 10

        # Set scan parameters
        arr = inparams
        if type(inparams) is not list:
            arr = [inparams]
            
        for inparam in arr:
            if inparam == 'hottail':
                self.addScanParameter(name='hottailgrid.pgrid.np',  f=_CS_setiHottailNp,  baselineValue=settings.hottailgrid.pgrid.np, scanUntilConvergence=scanUntilConvergence)
                self.addScanParameter(name='hottailgrid.xigrid.nxi', f=_CS_setiHottailNxi, baselineValue=settings.hottailgrid.xigrid.nxi, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'runaway':
                self.addScanParameter(name='runawaygrid.pgrid.np',  f=_CS_setiRunawayNp,  baselineValue=settings.runawaygrid.pgrid.np, scanUntilConvergence=scanUntilConvergence)
                self.addScanParameter(name='runawaygrid.xigrid.nxi', f=_CS_setiRunawayNxi, baselineValue=settings.runawaygrid.xigrid.nxi, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'hottailgrid.np':
                self.addScanParameter(name='hottailgrid.pgrid.np',  f=_CS_setiHottailNp,  baselineValue=settings.hottailgrid.pgrid.np, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'hottailgrid.nxi':
                self.addScanParameter(name='hottailgrid.xigrid.nxi', f=_CS_setiHottailNxi, baselineValue=settings.hottailgrid.xigrid.nxi, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'runawaygrid.np':
                self.addScanParameter(name='runawaygrid.pgrid.np',  f=_CS_setiRunawayNp,  baselineValue=settings.runawaygrid.pgrid.np, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'runawaygrid.nxi':
                self.addScanParameter(name='runawaygrid.xigrid.nxi', f=_CS_setiRunawayNxi, baselineValue=settings.runawaygrid.xigrid.nxi, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'nr':
                self.addScanParameter(name='radialgrid.nr', f=_CS_setiNr, baselineValue=settings.radialgrid.nr, scanUntilConvergence=scanUntilConvergence)
            elif inparam == 'nt':
                self.addScanParameter(name='timestep.nt', f=_CS_setiNt, baselineValue=settings.timestep.nt, scanUntilConvergence=scanUntilConvergence)
            elif inparam is None: continue
            else:
                self.addScanParameter(name=inparam, scanUntilConvergence=scanUntilConvergence)
                #raise ConvergenceScanException("Unrecognized scan parameter set: '{}'.".format(inparam))

        # Set output parameters
        arr = outparams
        if type(outparams) is not list:
            arr = [outparams]

        for outparam in arr:
            if outparam is None: continue
            else:
                self.addOutputParameter(name=outparam)


    def addOutputParameter(self, name: str, f=None, reltol=1e-2):
        """
        Adds an output parameter to check convergence for.

        :param str name:     Name of output parameter (used as an identifier, but does not have to correspond to the parameter's actual name in DREAM).
        :param function f:   A function which, given a DREAMOutput object, returns a single float value corresponding to this output parameter.
        :param float reltol: Relative tolerance to demand if 'scanUntilConvergence' is ``True`` for any of the scan parameters.
        """
        if f is None:
            # Try to get parameter by name
            f = lambda do : _CS_getoByName(do, name)

        self.outputParameters[name] = {'f': f, 'reltol': reltol}


    def addScanParameter(self, name: str, f=None, baselineValue=None, scanUntilConvergence=False, nvalues=3, startindex=-1):
        """
        Adds an input parameter to scan in.

        :param str name:                  Name of parameter (used as an identifier, but does not have to correspond to the parameter's actual name in DREAM).
        :param function f:                A function which, given an index, a ``DREAMSettings`` object and a baseline value, updates the scan parameter in the settings object. The index can take both positive and negative values, with ``0`` corresponding to the baseline value (negative values thus correspond to *lower resolution* while positive values correspond to *higher resolution*). The function should return a tuple consisting of the modified settings object (which may be the same as the input object) and the value representing the changes made to the ``DREAMSettings`` object (for identification purposes in plots).
        :param baselineValue:             Baseline value of the parameter (for passing on to ``f``).
        :param bool scanUntilConvergence: If ``True``, does not limit the number of runs to do and updates the value of each parameter until the output parameter changes less than the given tolerance.
        :param int nvalues:               Number of values to scan over. Ignored if ``scanUntilConvergence = True``.
        :param int startindex:            First index to run from. Default: ``-1`` (so that runs are ``-1``, ``0``, ``1``, ..., ``nvalues-2``)
        """
        if f is None:
            f = lambda idx, ds, v : _CS_setiByName(idx, ds, v, name)
        if baselineValue is None:
            baselineValue = self._getBaselineValue(name)

        self.scanParameters[name] = {'f': f, 'baseline': baselineValue, 'scanUntilConvergence': scanUntilConvergence, 'nvalues': int(nvalues), 'startindex': int(startindex)}
        self.result[name] = {}


    def _getBaselineValue(self, name):
        """
        Returns the baseline value for the named input parameter.
        """
        obj = _CS_getObjectByName(self.settings, name)

        if np.isscalar(obj):
            return obj
        else:
            raise ConvergenceScanException("Unrecognized type of input parameter '{}': {}.".format(name, type(obj)))


    def getOutputParameters(self):
        """
        Get a dictionary containing details about the output parameters used in the scan.

        :return: A dict which specifies the settings for the output parameters to consider as measures of convergence.
        :rtype: dict
        """
        oparams = {}
        for opname, op in self.outputParameters.items():
            oparams[opname] = {'reltol': op['reltol']}

        return oparams


    def _processOutput(self, index, scanParameter, scanValue, output):
        """
        Process the output of a single simulation.

        index:         Index of run.
        scanParameter: Name of scan parameter settings specifying the scan.
        output:        DREAMOutput object resulting from the simulation.

        RETURNS False if the relative variation in any output parameter,
        compared to the previous (or next) index, exceeds its relative
        tolerance, or if there is no previous or next scan index.
        Otherwise, returns True (meaning that the parameter is converged).
        """
        sp = self.scanParameters[scanParameter]
        converged = True

        if scanParameter not in self.result:
            self.result[scanParameter] = {}

        for opname, op in self.outputParameters.items():
            f = op['f']
            oval = f(output)

            # Store output value in result
            if opname not in self.result[scanParameter]:
                self.result[scanParameter][opname] = {'index': [], 'scanval': [], 'outval': []}

            self.result[scanParameter][opname]['index'].append(index)
            self.result[scanParameter][opname]['scanval'].append(scanValue)
            self.result[scanParameter][opname]['outval'].append(oval)

            # Check if simulation is converged
            if len(self.result[scanParameter][opname]['outval']) == 1:
                converged = False
            else:
                v = self.result[scanParameter][opname]['outval'][-2:]
                reltol = op['reltol']

                # We should really have an absolute tolerance for this...
                cv = False
                Delta = 1
                if v[1] == 0:
                    Delta = np.abs(v[0]) 
                    cv = (Delta < reltol)
                else:
                    Delta = np.abs(v[0]/v[1] - 1)
                    cv = (Delta < reltol)

                converged = converged and cv
                self._status('Output parameter {} {} converged in scan parameter {} (Delta = {})'.format(opname, 'is' if cv else 'is not', scanParameter, Delta))

        return converged


    def run(self):
        """
        Run the convergence scan.
        """
        self.result = {}

        # Run baseline case
        self._status(':: Running baseline case...')
        self.baselineOutput = runiface.runiface(self.settings, quiet=not self.verbose)

        # Iterate over scan parameters
        for scanParameter, sp in self.scanParameters.items():
            s = sp['startindex']

            if sp['scanUntilConvergence']:
                n = self.NMAX + s
                i = s
                converged = False

                while i < n and not converged:
                    output, scanValue = self._runScan(i, scanParameter)
                    converged = self._processOutput(i, scanParameter, scanValue, output)
                    i += 1
            else:
                n = sp['nvalues']
                for i in range(s, n+s):
                    output, scanValue = self._runScan(i, scanParameter)
                    self._processOutput(i, scanParameter, scanValue, output)
                    

    def _runScan(self, index, scanParameter):
        """
        Run an individual DREAM simulation corresponding to index
        'index' in the scan parameter 'scanParameter'.

        index:         Index in scan of this simulation.
        scanParameter: Name of scan parameter settings specifying the scan.
        """
        sp = self.scanParameters[scanParameter]
        self._status(':: Scan {} ({}/{}) in parameter {}'.format(index, index-sp['startindex']+1, sp['nvalues'], scanParameter))

        # Skip the baseline case
        if index == 0:
            self._status(':: Skipping baseline case')
            return self.baselineOutput, sp['baseline']

        f = sp['f']
        # Copy DREAMSettings object
        ns = DREAMSettings(self.settings, chain=False)

        # Modify the settings
        ns, scanValue = f(index, ns, sp['baseline'])

        return runiface.runiface(ns), scanValue


    def save(self, filename):
        """
        Saves this convergence scan to an HDF5 file.

        :param str filename: Name of file to save scan results to.
        """
        DREAMIO.SaveDictAsHDF5(filename=filename, data=self.todict())


    def setVerbose(self, verbose=True):
        """
        If verbose is ``True``, the scan will print progress messages to stdout when running.

        :param bool verbose: Value to set verbosity to.
        """
        self.verbose = verbose


    def _status(self, msg):
        if self.verbose:
            print(msg)


    def todict(self):
        """
        Converts the results of this scan to a dictionary object which can easily be saved to file.
        """
        oparams = self.getOutputParameters()
        
        return {
            'result': self.result,
            'outputParameters': oparams
        }


class ConvergenceScanException(Exception):
    def __init__(self, msg):
        super(Exception, self).__init__(msg)


def _CS_getObjectByName(do, name, paramType='input'):
    lst = []
    if '.' in name: lst = name.split('.')
    else: lst = [name]

    obj = do
    for o in lst:
        if o in obj:
            obj = obj[o]
        else:
            raise ConvergenceScanException("Unrecognized {} parameter specified: '{}'.".format(paramType, o))

    return obj


def _CS_setObjectByName(ds, name, val, paramType='input'):
    lst = []
    if '.' in name: lst = name.split('.')
    else: lst = [name]

    obj = ds
    for o in lst[:-1]:
        if o in obj:
            obj = obj[o]
        else:
            raise ConvergenceScanException("Unrecognized {} parameter specified: '{}'.".format(paramType, name))

    if not np.isscalar(obj[lst[-1]]):
        raise ConvergenceScanException("The input parameter '{}' is of an unrecognized type: {}.".format(name, type(obj)))

    obj.__dict__[lst[-1]] = val


# Helper functions for updating resolution parameters
def _CS_setiByName(index: int, settings: DREAMSettings, baseline, name: str):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    _CS_setObjectByName(settings, name, val, 'input')
    return settings, val

def _CS_setiHottailNp(index: int, settings: DREAMSettings, baseline):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    settings.hottailgrid.setNp(val)
    return settings, val

def _CS_setiHottailNxi(index: int, settings: DREAMSettings, baseline):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    settings.hottailgrid.setNxi(val)
    return settings, val

def _CS_setiRunawayNp(index: int, settings: DREAMSettings, baseline):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    settings.runawaygrid.setNp(val)
    return settings, val

def _CS_setiRunawayNxi(index: int, settings: DREAMSettings, baseline):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    settings.runawaygrid.setNxi(val)
    return settings, val

def _CS_setiNr(index: int, settings: DREAMSettings, baseline):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    settings.radialgrid.setNr(val)
    return settings, val

def _CS_setiNt(index: int, settings: DREAMSettings, baseline):
    val = max(1,int(np.round(baseline * np.float_power(2, index))))
    settings.timestep.setNt(val)
    return settings, val


def _CS_getoByName(do: DREAMOutput, param: str) -> float:
    obj = _CS_getObjectByName(do, param, 'output')

    # Get last element of array (this makes most
    # sense for time-evolving scalar data, as this
    # selection then corresponds to the final value
    # of the scalar parameter)
    arr = obj[:]
    slc = tuple([-1] * arr.ndim)

    return arr[slc]

