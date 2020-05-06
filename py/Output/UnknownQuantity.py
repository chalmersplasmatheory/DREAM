

import numpy as np
import matplotlib.pyplot as plt


class UnknownQuantity:
    

    def __init__(self, name, data, grid):
        """
        Constructor.

        name: Name of unknown.
        data: Data of unknown.
        grid: Grid used for the DREAM simulation.
        """
        self.name = name
        self.data = data
        self.grid = grid


    def __getitem__(self, key):
        """
        Direct access to 'data' dict.
        """
        return self.data[key]


    def getName(self): return self.name


    def getData(self): return self.data
    

    def getTeXName(self):
        return self.name.replace('_', r'\_')


