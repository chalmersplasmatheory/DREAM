
import numpy as np
import matplotlib.pyplot as plt


class OtherQuantity:
    

    def __init__(self, name, data, description, grid, output, momentumgrid=None):
        """
        Constructor.

        name:         Name of quantity.
        data:         Data of quantity (raw, as dict from output file).
        description:  String describing this quantity.
        grid:         Grid that was used for the DREAM simulation.
        output:       Parent DREAMOutput object.
        momentumgrid: Momentum grid associated with quantity.
        """
        self.name         = name
        self.data         = data
        self.description  = description
        self.grid         = grid
        self.output       = output
        self.momentumgrid = momentumgrid


    def __getitem__(self, key):
        """
        Direct access to 'data' dict.
        """
        return self.data[key]


    def __repr__(self):
        return self.__str__()


    def __str__(self):
        return "({}) {}".format(self.name, self.description)


    def getName(self): return self.name


    def getData(self): return self.data


    def getTeXName(self):
        return self.name.replace('_', r'\_')


    def getTeXIntegralName(self):
        return 'Integrated '+self.getTeXName()


