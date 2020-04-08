# Equation system settings object
###################################

import numpy as np

from Equations.ColdElectrons import ColdElectrons


class EquationSystem:
    
    def __init__(self):
        """
        Constructor.
        """
        self.unknowns = {
            'n_cold': ColdElectrons()
        }


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this EquationSystem object.
        """
        data = {}

        for u, obj in self.unknowns.items():
            data[u] = obj.todict()

        return data


    def verifySettings(self):
        """
        Verify that all unknowns have been properly configured
        and that all settings are consistent.
        """
        for u in self.unknowns:
            u.verifySettings()

