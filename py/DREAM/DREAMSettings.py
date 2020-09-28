#
# An object representing the settings passed when running DREAM.
# ###############################################################

import copy
import numpy as np
import os
from . import DREAMIO as DREAMIO

# Settings objects
from .Settings.CollisionHandler import CollisionHandler
from .Settings.EquationSystem import EquationSystem
from .Settings.MomentumGrid import MomentumGrid
from .Settings.OtherQuantities import OtherQuantities
from .Settings.Output import Output
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
        self.init = {}

        self.addSetting('collisions', CollisionHandler())
        self.addSetting('hottailgrid', MomentumGrid('hottailgrid'))
        self.addSetting('other', OtherQuantities())
        self.addSetting('output', Output())
        self.addSetting('radialgrid', RadialGrid())
        self.addSetting('runawaygrid', MomentumGrid('runawaygrid'))
        self.addSetting('solver', Solver())
        self.addSetting('timestep', TimeStepper())

        # Should be defined last as it may need access to the
        # objects created above...
        self.addSetting('eqsys', EquationSystem(settings=self))

        if filename is not None:
            if type(filename) == str:
                self.load(filename, path=path)
            elif type(filename) == DREAMSettings:
                self.fromdict(filename.todict())

                self.fromOutput(filename.output.filename)
                self.output.setFilename('output.h5')

    
    def __contains__(self, item):
        """
        Overriding the Python 'in' keyword.
        """
        return (item in self.settings)


    def __getitem__(self, index):
        """
        Retrieves a parameter by name.
        """
        return self.settings[index]


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


    def fromdict(self, data, filename='<dictionary>'):
        sets  = list(self.settings.keys())
        other = ['init']

        for key in data:
            # Warn about unrecognized settings
            if key in sets:
                # Remove from list of not-found settings
                sets.remove(key)
                # Set settings
                if type(self.settings[key]) == MomentumGrid:
                    self.settings[key].fromdict(key, data[key])
                else:
                    self.settings[key].fromdict(data[key])
            elif key in other:
                # Remove from list of not found
                other.remove(key)

                # Set settings
                setattr(self, key, data[key])
            else:
                print("WARNING: Unrecognized setting '{}' found in '{}'.".format(key, filename))
                continue

        # Warn about missing settings
        missing = sets+other
        if len(missing) > 0:
            for s in missing:
                print("WARNING: Setting '{}' not specified in '{}'.".format(s, filename))


    def fromOutput(self, filename, relpath=False, ignore=list(), timeindex=-1):
        """
        Specify that the simulation should be initialized from the
        DREAM output stored in the named file. Some unknown quantities
        can be ignored and initialized conventionally by adding them
        to the 'ignore' list.

        filename:  Name of file to load output from.
        relpath:   If 'True', forces the path to 'filename' to be encoded
                   as a relative path. Otherwise, the absolute path to the
                   named file will be calculated and given to the settings
                   object. (default: False)
        ignore:    List of unknown quantities to initialize as usual, and
                   thus NOT from the specified output file.
        timeindex: Index of time point to use for initializing the simulation.
        """
        # Input as relative or absolute path?
        if relpath:
            fname = filename
        else:
            fname = os.path.abspath(filename)

        if type(ignore) == str:
            ignore = [ignore]
        elif type(ignore) != list:
            raise DREAMException("Unrecognized type of argument 'ignore'. Expected list of strings.")

        self.init['fromfile']   = fname
        self.init['ignorelist'] = ignore
        self.init['timeindex']  = timeindex


    def load(self, filename, path=""):
        """
        Load a DREAMSettings object from the named HDF5 file.
        'path' specifies the path within the HDF5 file where
        the DREAMSettings object is stored.
        """
        data = DREAMIO.LoadHDF5AsDict(filename, path=path)
        self.fromdict(data, filename=filename)


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

        data['init'] = {}

        if ('ignorelist' in self.init) and (len(self.init['ignorelist']) > 0):
            data['init']['eqsysignore'] = ';'.join(self.init['ignorelist'])

        if 'timeindex' in self.init:
            data['init']['filetimeindex'] = self.init['timeindex']
        if 'fromfile' in self.init:
            data['init']['fromfile'] = self.init['fromfile']

        return data


    def verifySettings(self):
        """
        Verify that the DREAM run has been correctly configured
        and that all settings are consistent.
        """
        for _, setting in self.settings.items():
            setting.verifySettings()

        if ('fromfile' in self.init) and (self.init['fromfile'] != ''):
            # Verify that the file exists...
            if not os.path.exists(self.init['fromfile']):
                print("WARNING: The output file from which to initialize '{}' does not exist.".format(self.init['fromfile']))
        

