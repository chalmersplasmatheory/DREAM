# Equation system settings object
###################################

import numpy as np
from .. DREAMException import DREAMException

from .Equations.ColdElectrons import ColdElectrons
from .Equations.ColdElectronTemperature import ColdElectronTemperature
from .Equations.ElectricField import ElectricField
from .Equations.HotElectronDistribution import HotElectronDistribution
from .Equations.Ions import Ions
from .Equations.EquationException import EquationException


class EquationSystem:
    
    def __init__(self):
        """
        Constructor.
        """
        self.unknowns = {}
        self.addUnknown('E_field', ElectricField())
        self.addUnknown('f_hot', HotElectronDistribution())
        self.addUnknown('n_cold', ColdElectrons())
        self.addUnknown('n_i', Ions())
        self.addUnknown('T_cold', ColdElectronTemperature())


    def addUnknown(self, name, obj):
        """
        Add an unknown to this object. This adds the unknown to
        the 'unknowns' list, in addition to making it accessible
        through usual "dot" notation (i.e. you can access it either
        as "self.myUnknown" or "self.unknowns['myUnknown']")

        name: Name of unknown object to add.
        obj:  Unknown object to add.
        """
        setattr(self, name, obj)
        self.unknowns[name] = obj


    def fromdict(self, data):
        """
        Sets the options of this object from a dictionary.
        """
        sets = list(self.unknowns.keys())

        for key in data:
            # Warn about unrecognized settings
            if key not in sets:
                print("WARNING: Unrecognized setting '{}'.".format(key))
                continue

            # Remove from list of not-found settings
            sets.remove(key)
            # Set settings
            try:
                self.unknowns[key].fromdict(data[key])
            except EquationException as ex:
                print("WARNING: {}".format(ex))

        # Warn about missing settings
        if len(sets) > 0:
            for s in sets:
                print("WARNING: Settings for unknown '{}' not specified.".format(s))


    def todict(self, verify=True):
        """
        Returns a Python dictionary containing all settings of
        this EquationSystem object.
        """
        if verify:
            self.verifySettings()

        data = {}

        for u, obj in self.unknowns.items():
            data[u] = obj.todict()

        return data


    def verifySettings(self):
        """
        Verify that all unknowns have been properly configured
        and that all settings are consistent.
        """
        for _, u in self.unknowns.items():
            u.verifySettings()

