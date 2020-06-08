#

import numpy as np
from . FluidQuantity import FluidQuantity
from . HotElectronDistributionFunction import HotElectronDistributionFunction
from . IonHandler import IonHandler
from . UnknownQuantity import UnknownQuantity


class EquationSystem:


    SPECIAL_TREATMENT = {
        # List of unknown quantities with their own classes
        'E_field': FluidQuantity,
        'f_hot':   HotElectronDistributionFunction,
        'j_hot':   FluidQuantity,
        'n_cold':  FluidQuantity,
        'n_hot':   FluidQuantity,
        'n_i':     IonHandler,
        'n_re':    FluidQuantity,
        'n_tot':   FluidQuantity,
        'T_cold':  FluidQuantity
    }
    

    def __init__(self, unknowns=None, grid=None, output=None):
        """
        Constructor.

        unknowns: List of unknowns in the equation system (with data).
        """
        self.grid = grid
        self.unknowns = {}
        self.output = output

        if unknowns is not None:
            self.setUnknowns(unknowns)


    def __getitem__(self, index):
        """
        Direct access by name to the list of unknowns.
        """
        return self.unknowns[index]


    def keys(self): return self.unknowns.keys()


    def getUnknownNames(self):
        """
        Get a list with the names of all unknowns in
        the equation system.
        """
        return list(self.unknowns.keys())


    def setGrid(self, grid):
        """
        Sets the grid used for the DREAM simulation.
        """
        self.grid = grid


    def setUnknown(self, name, data):
        """
        Add the given unknown to this equation system.

        name: Name of the unknown.
        data: Data for the unknown (raw, as a dict from the output file).
        """
        if name in self.SPECIAL_TREATMENT:
            o = self.SPECIAL_TREATMENT[name](name=name, data=data, grid=self.grid, output=self.output)
        else:
            o = UnknownQuantity(name=name, data=data, grid=self.grid, output=self.output)

        setattr(self, name, o)
        self.unknowns[name] = o
        

    def setUnknowns(self, unknowns):
        """
        Add a list of unknowns to this equation system.
        """
        for uqn in unknowns:
            self.setUnknown(name=uqn, data=unknowns[uqn])


