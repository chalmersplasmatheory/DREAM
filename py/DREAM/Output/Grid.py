
import numpy as np


class Grid:
    

    def __init__(self, grid=None):
        """
        Constructor.

        grid: Grid data from DREAM output.
        """
        self.t = None
        self.r = None
        self.hottail = None
        self.runaway = None

        if grid is not None:
            self.setGrid(grid)


    def setGrid(self, grid):
        """
        Set grid data based on output from DREAM.
        """
        self.t = grid['t']
        self.r = grid['r']

        if 'hottail' in grid:
            self.hottail = grid['hottail']

        if 'runaway' in grid:
            self.runaway = grid['runaway']

