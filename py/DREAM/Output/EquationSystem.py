#

import numpy as np
from . CurrentDensity import CurrentDensity
from . ElectricField import ElectricField
from . FluidQuantity import FluidQuantity
from . HotElectronDistributionFunction import HotElectronDistributionFunction
from . RunawayElectronDistributionFunction import RunawayElectronDistributionFunction
from . IonHandler import IonHandler
from . ScalarQuantity import ScalarQuantity
from . UnknownQuantity import UnknownQuantity


class EquationSystem:


    SPECIAL_TREATMENT = {
        # List of unknown quantities with their own classes
        'E_field':      ElectricField,
        'f_hot':        HotElectronDistributionFunction,
        'f_re':         RunawayElectronDistributionFunction,
        'I_p':          ScalarQuantity,
        'I_wall':       ScalarQuantity,
        'j_hot':        CurrentDensity,
        'j_ohm':        CurrentDensity,
        'j_re':         CurrentDensity,
        'j_tot':        CurrentDensity,
        'n_cold':       FluidQuantity,
        'n_hot':        FluidQuantity,
        'n_i':          IonHandler,
        'n_re':         FluidQuantity,
        'n_tot':        FluidQuantity,
        'psi_edge':     ScalarQuantity,
        'psi_p':        FluidQuantity,
        'psi_wall':     ScalarQuantity,
        'S_particle':   FluidQuantity,
        'tau_coll':     FluidQuantity,
        'T_cold':       FluidQuantity,
        'V_loop_w':     ScalarQuantity,
        'W_cold':       FluidQuantity
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


    def __contains__(self, item):
        """
        Overrides the Python 'in' operator.
        """
        return (item in self.unknowns)


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


    def setUnknown(self, name, data, attr):
        """
        Add the given unknown to this equation system.

        name: Name of the unknown.
        data: Data for the unknown (raw, as a dict from the output file).
        attr: List of attributes set to this unknown in the output file.
        """
        if name in self.SPECIAL_TREATMENT:
            o = self.SPECIAL_TREATMENT[name](name=name, data=data, attr=attr, grid=self.grid, output=self.output)
        else:
            o = UnknownQuantity(name=name, data=data, attr=attr, grid=self.grid, output=self.output)

        setattr(self, name, o)
        self.unknowns[name] = o
        

    def setUnknowns(self, unknowns):
        """
        Add a list of unknowns to this equation system.
        """
        for uqn in unknowns:
            # Skip attributes
            if uqn[-2:] == '@@': continue

            attr = []
            if uqn+'@@' in unknowns:
                attr = unknowns[uqn+'@@']

            self.setUnknown(name=uqn, data=unknowns[uqn], attr=attr)


