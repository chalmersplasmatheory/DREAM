#

import numpy as np
from Output.FluidQuantity import FluidQuantity
from Output.UnknownQuantity import UnknownQuantity


class EquationSystem:


    SPECIAL_TREATMENT = {
        # List of unknown quantities with their own classes
        'E_field': FluidQuantity,
        'T_cold':  FluidQuantity,
        'n_cold':  FluidQuantity,
        'n_hot':   FluidQuantity,
        'n_re':    FluidQuantity,
        'n_tot':   FluidQuantity
    }
    

    def __init__(self, unknowns=None, grid=None):
        """
        Constructor.

        unknowns: List of unknowns in the equation system (with data).
        """
        self.grid = grid
        self.unknowns = {}

        if unknowns is not None:
            self.setUnknowns(unknowns)


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
        data: Data for the unknown.
        """
        if name in self.SPECIAL_TREATMENT:
            o = self.SPECIAL_TREATMENT[name](name=name, data=data, grid=self.grid)
        else:
            o = UnknownQuantity(name=name, data=data, grid=self.grid)

        setattr(self, name, o)
        self.unknowns[name] = o
        

    def setUnknowns(self, unknowns):
        """
        Add a list of unknowns to this equation system.
        """
        for uqn in unknowns:
            self.setUnknown(name=uqn, data=unknowns[uqn])


