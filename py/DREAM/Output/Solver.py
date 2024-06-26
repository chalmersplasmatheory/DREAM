# Output handler for solver statistics

import numpy as np
from ..DREAMException import DREAMException
from ..DataObject import DataObject

import DREAM.Settings.Solver


class Solver:
    

    def __init__(self, solverdata=None, output=None):
        """
        Constructor.
        """
        self.solverdata = None
        self.output = output
        self.warnings = []

        if 'warnings' in solverdata:
            self.warnings = solverdata['warnings'][:].split(';')[:-1]


    def __contains__(self, item):
        """
        Overrides the Python 'in' operator.
        """
        return (item in self.__dict__)


    def __getitem__(self, index):
        """
        Direct access by name to the timing information.
        """
        return self.__dict__[index]


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def plot(self):
        """
        Generic plotting routine
        """
        print("Method 'plot()' not implemented for this solver.")


