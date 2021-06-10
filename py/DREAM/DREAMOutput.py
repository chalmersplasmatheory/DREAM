#
# This object represents the output of a DREAM simulation.
# ###########################################################

import copy
import numpy as np
import os
import DREAM.DREAMIO as DREAMIO

from .DREAMSettings import DREAMSettings
from .Output.EquationSystem import EquationSystem
from .Output.Grid import Grid
from .Output.IonMetaData import IonMetaData
from .Output.OtherQuantityHandler import OtherQuantityHandler
from .Output.SolverLinear import SolverLinear
from .Output.SolverNonLinear import SolverNonLinear
from .Output.Timings import Timings

from .Settings import Solver as SettingsSolver


class DREAMOutput:
    

    def __init__(self, filename=None, path="", lazy=True):
        """
        Construct a new ``DREAMOutput`` object. If ``filename`` is given, the
        object is read from the (HDF5) file with that name. If ``path`` is also
        given, this is used to locate the group in the file which contains the
        output.

        :param str filename: Name of file to load output from.
        :param str path:     Path to group in HDF5 file containing the output.
        """

        # Default
        self.eqsys = None
        self.grid = None
        self.ionmeta = None
        self.other = None
        self.settings = None
        self.solver = None
        self.timings = None

        self.filename = None
        self.filesize = 0
        self.h5handle = None

        if filename is not None:
            self.load(filename=filename, path=path, lazy=lazy)


    def __contains__(self, item):
        """
        Overriding the Python ``in`` operator and allows to check for the
        existence of properties of this object.
        """
        return (item in self.__dict__)


    def __getitem__(self, index):
        """
        Retrieves a parameter by name.
        """
        return self.__dict__[index]


    def close(self):
        """
        Close the associated HDF5 File object.
        """
        if self.h5handle is not None:
            self.h5handle.close()


    def _getsolver(self, solverdata, output):
        if 'type' in solverdata:
            if solverdata['type'] == SettingsSolver.LINEAR_IMPLICIT:
                return SolverLinear(solverdata, output)
            elif solverdata['type'] == SettingsSolver.NONLINEAR:
                return SolverNonLinear(solverdata, output)
        else:
            print('WARNING: Invalid solver data given.')
            return None


    def load(self, filename, path="", lazy=True):
        """
        Loads DREAM output from the specified file. If 'path' is
        given, this indicates which group path in the file to load
        the output from.

        :param str filename: Name of file to load output from.
        :param str path:     Path to subsect of HDF5 file containing DREAM output.
        :param bool lazy:    If ``True``, allows the file to be read lazily (on-demand) by return h5py DataSet objects instead of the actual data (wrapped in a DREAM.DataObject).  This can greatly reduce load times, but may complicate typing slightly. Note also that the HDF5 file will be locked for as long as the Python interpreter is running.
        """
        self.filename = filename

        od, self.h5handle, self.filesize = DREAMIO.LoadHDF5AsDict(filename, path=path, returnhandle=True, returnsize=True, lazy=lazy)

        if 'grid' in od:
            self.grid = Grid(od['grid'])
        else:
            print("WARNING: No grid found in '{}'.".format(filename))
        
        if 'ionmeta' in od:
            self.ionmeta = IonMetaData(od['ionmeta'])
        else:
            print("WARNING: No ion meta data found in '{}'.".format(filename))

        # Equation system should be loaded last, because it
        # may depend on previously loaded sections
        if 'eqsys' in od:
            self.eqsys = EquationSystem(od['eqsys'], grid=self.grid, output=self)
        else:
            print("WARNING: No equation system found in '{}'.".format(filename))
            
        
        # Load "other" quantities (i.e. quantities which are not part of
        # the equation system, but may still be interesting to know the
        # evolution of; this include collision frequencies, bounce averages
        # and more)
        if 'other' in od:
            self.other = OtherQuantityHandler(od['other'], grid=self.grid, output=self)

        # Load settings for the run
        if 'settings' in od:
            self.settings = DREAMSettings(settings=od['settings'])

        # Solver statistics
        if 'solver' in od:
            self.solver = self._getsolver(od['solver'], output=self)

        # Timing information
        if 'timings' in od:
            self.timings = Timings(od['timings'], output=self)


    def getFileSize(self):
        """
        Returns the size in bytes of the output file.
        """
        return self.filesize


    def getFileSize_s(self, includeBytes=False):
        """
        Returns the size in bytes of the output file as a
        nicely formatted string.
        """
        # Ensure functionality of future super-sized DREAM runs...
        unit = ['B', 'kiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB']

        i = 0
        fs = self.filesize
        while fs > 1024:
            fs /= 1024
            i += 1

        s = '{:.2f} {}'.format(fs, unit[i])

        if includeBytes and i > 0:
            s += ' ({} bytes)'.format(self.filesize)

        return s


