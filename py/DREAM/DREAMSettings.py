#
# An object representing the settings passed when running DREAM.
# ###############################################################

import copy
import numpy as np
from . import DREAMIO as DREAMIO

# Settings objects
from .Settings.EquationSystem import EquationSystem
from .Settings.MomentumGrid import MomentumGrid
from .Settings.OtherQuantities import OtherQuantities
from .Settings.RadialGrid import RadialGrid
from .Settings.Solver import Solver
from .Settings.TimeStepper import TimeStepper


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

        self.addSetting('eqsys', EquationSystem())
        self.addSetting('hottailgrid', MomentumGrid('hottailgrid'))
        self.addSetting('other', OtherQuantities())
        self.addSetting('radialgrid', RadialGrid())
        self.addSetting('runawaygrid', MomentumGrid('runawaygrid'))
        self.addSetting('solver', Solver())
        self.addSetting('timestep', TimeStepper())

        if filename is not None:
            self.load(filename, path=path)

    
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


    def load(self, filename, path=""):
        """
        Load a DREAMSettings object from the named HDF5 file.
        'path' specifies the path within the HDF5 file where
        the DREAMSettings object is stored.
        """
        data = DREAMIO.LoadHDF5AsDict(filename, path=path)
        sets = list(self.settings.keys())

        for key in data:
            # Warn about unrecognized settings
            if key not in sets:
                print("WARNING: Unrecognized setting '{}' found in '{}'.".format(key, filename))
                continue

            # Remove from list of not-found settings
            sets.remove(key)
            # Set settings
            if type(self.settings[key]) == MomentumGrid:
                self.settings[key].fromdict(key, data[key])
            else:
                self.settings[key].fromdict(data[key])
                
        # Warn about missing settings
        if len(sets) > 0:
            for s in sets:
                print("WARNING: Setting '{}' not specified in '{}'.".format(s, filename))


    def save(self, filename):
        """
        Save this settings object to the specified file.

        filename: Name of file to save settings to (the file will
                  be overwritten if it exists).
        """
        DREAMIO.SaveDictAsHDF5(filename, self.todict())


    def todict(self, verify=True):
        """
        Returns the settings in this object as a Python dictionary.
        """
        if verify:
            self.verifySettings()

        data = {}
        for key, setting in self.settings.items():
            data[key] = setting.todict(verify=False)

        return data


    def verifySettings(self):
        """
        Verify that the DREAM run has been correctly configured
        and that all settings are consistent.
        """
        for _, setting in self.settings.items():
            setting.verifySettings()
        

