# Equation system settings object
###################################

import copy
import numpy as np
from .. DREAMException import DREAMException

from .Equations.ColdElectrons import ColdElectrons
from .Equations.ColdElectronTemperature import ColdElectronTemperature
from .Equations.ElectricField import ElectricField
from .Equations.HotElectronDistribution import HotElectronDistribution
from .Equations.Ions import Ions
from .Equations.RunawayElectrons import RunawayElectrons
from .Equations.EquationException import EquationException


# List of names of unknown quantities in DREAM. This list can be
# used across the interface to validate names of unknowns.
UNKNOWNS = [
    'E_field', 'f_hot', 'f_re', 'n_i', 'I_p', 'I_wall',
    'j_hot', 'j_ohm', 'j_re', 'j_tot', 'n_cold', 'n_hot',
    'n_re', 'n_tot', 'psi_p', 'psi_wall', 'psi_edge',
    'T_cold', 'V_loop_w', 'W_cold'
]


class EquationSystem:
    
    def __init__(self, settings):
        """
        Constructor.

        settings: Parent settings object.
        """
        # Parent settings object
        self.settings = settings

        self.unknowns = list()
        self.addUnknown('E_field', ElectricField(settings=settings))
        self.addUnknown('f_hot', HotElectronDistribution(settings=settings))
        self.addUnknown('n_cold', ColdElectrons(settings=settings))
        self.addUnknown('n_i', Ions(settings=settings))
        self.addUnknown('n_re', RunawayElectrons(settings=settings))
        self.addUnknown('T_cold', ColdElectronTemperature(settings=settings))


    def __getitem__(self, name):
        """
        Get UnknownQuantity by name.
        """
        return getattr(self, name)


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
        self.unknowns.append(name)


    def fromdict(self, data):
        """
        Sets the options of this object from a dictionary.
        """
        sets = copy.copy(self.unknowns)

        for key in data:
            # Warn about unrecognized settings
            if key not in sets:
                print("WARNING: Unrecognized setting '{}'.".format(key))
                continue

            # Remove from list of not-found settings
            sets.remove(key)
            # Set settings
            try:
                self[key].fromdict(data[key])
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

        for u in self.unknowns:
            data[u] = self[u].todict()

        return data


    def verifySettings(self):
        """
        Verify that all unknowns have been properly configured
        and that all settings are consistent.
        """
        for u in self.unknowns:
            self[u].verifySettings()


