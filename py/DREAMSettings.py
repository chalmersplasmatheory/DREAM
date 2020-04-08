#
# An object representing the settings passed when running DREAM.
# ###############################################################

import copy
import numpy as np
import DREAMIO

# Settings objects
from TimeStepper import TimeStepper


class DREAMSettings:
    
    TIMESTEP_TYPE_CONSTANT = 1
    
    def __init__(self, filename=None, path=""):
        """
        Construct a new DREAMSettings object. If 'filename' is given,
        the object is read from the (HDF5) file with that name.
        If 'path' is also given, this is used to locate the group
        in the file which contains the settings. 

        filename: Name of the file to load settings from.
        path:     Path to group in HDF5 file containing the settings.
        """

        # Defaults
        self.settings = {}

        self.addSetting('equationsystem', EquationSystem())
        self.addSetting('timestep', TimeStepper())

    
    def addSetting(self, name, obj):
        """
        Add a setting to this object. This adds the setting to
        the 'settings' list, in addition to making it accessible
        through usual "dot" notation (i.e. you can access it either
        as "self.mySetting" or "self.settings['mySetting']")

        name: Name of settings object to add.
        obj:  Settings object to add.
        """
        setattr(self, name, obj)
        self.settings[name] = obj


    def save(self, filename):
        """
        Save this settings object to the specified file.

        filename: Name of file to save settings to (the file will
                  be overwritten if it exists).
        """
        DREAMIO.SaveDictAsHDF5(filename, self.todict())


    def todict(self):
        """
        Returns the settings in this object as a Python dictionary.
        """
        data = {}
        for key, setting in self.settings.items():
            data[key] = setting.todict()

        return data
        

