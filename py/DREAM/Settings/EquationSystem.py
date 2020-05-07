# Equation system settings object
###################################

import numpy as np
from DREAM.DREAMException import DREAMException

from DREAM.Settings.Equations.ColdElectrons import ColdElectrons
from DREAM.Settings.Equations.ColdElectronTemperature import ColdElectronTemperature
from DREAM.Settings.Equations.ElectricField import ElectricField
from DREAM.Settings.Equations.HotElectronDistribution import HotElectronDistribution
from DREAM.Settings.Equations.Ions import Ions


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

