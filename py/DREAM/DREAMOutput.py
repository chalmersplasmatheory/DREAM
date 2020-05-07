#
# This object represents the output of a DREAM simulation.
# ###########################################################

import copy
import numpy as np
import DREAM.DREAMIO as DREAMIO

from DREAM.DREAMSettings import DREAMSettings
from DREAM.Output.EquationSystem import EquationSystem
from DREAM.Output.Grid import Grid


class DREAMOutput:
    

    def __init__(self, filename=None, path=""):
        """
        Construct a new DREAMOutput object. If 'filename' is given,
        the object is read from the (HDF5) file with that name.
        If 'path' is also given, this is used to locate the group
        in the file which contains the settings.

        filename: Name of file to load output from.
        path:     Path to group in HDF5 file containing the output.
        """

        # Default
        self.eqsys = None
        self.grid = None
        self.settings = None

        if filename is not None:
            self.load(filename=filename, path=path)


    def load(self, filename, path=""):
        """
        Loads DREAM output from the specified file. If 'path' is
        given, this indicates which group path in the file to load
        the output from.

        filename: Name of file to load output from.
        path:     Path to output in HDF5 file.
        """
        os = DREAMIO.LoadHDF5AsDict(filename, path=path)

        if 'grid' in os:
            self.grid = Grid(os['grid'])
        else:
            print("WARNING: No grid found in '{}'.".format(filename))

        if 'eqsys' in os:
            self.eqsys = EquationSystem(os['eqsys'], grid=self.grid)
        else:
            print("WARNING: No equation system found in '{}'.".format(filename))

        if 'settings' in os:
            # TODO Make this work!
            self.settings = DREAMSettings(settings=os['settings'])


