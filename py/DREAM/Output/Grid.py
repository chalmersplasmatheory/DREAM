
import numpy as np
import DREAM.Settings.MomentumGrid as MomentumGrid


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


    def __str__(self):
        """
        Convert this object to a string.
        """
        s  = "DREAM Grid Object\n"
        s += "   t: {} elements from {} to {}\n".format(self.t.size, self.t[0], self.t[-1])
        s += "   r: {} elements from {} to {}\n".format(self.r.size, self.r[0], self.r[-1])
        s += "   hottail   ({})\n{}".format(self.strType(self.hottail) if self.hottail is not None else 'DISABLED', self.strGrid(self.hottail))
        s += "   runaway   ({})\n{}".format(self.strType(self.runaway) if self.runaway is not None else 'DISABLED', self.strGrid(self.runaway))

        return s


    def strType(self, grid):
        """
        Returns the grid type as a string representation
        """
        if grid is None: return ""

        if grid['type'] == MomentumGrid.MOMENTUMGRID_TYPE_PXI:
            return 'p/xi'
        elif grid['type'] == MomentumGrid.MOMENTUMGRID_TYPE_PPARPPERP:
            return 'ppar/pperp'
        else:
            return '<unknown-type>'


    def strGrid(self, grid):
        """
        Returns some information about the grid as a string.
        """
        if grid is None: return ""

        s = ""
        if grid['type'] == MomentumGrid.MOMENTUMGRID_TYPE_PXI:
            s += "      p:  {} elements from {} to {}\n".format(grid['p1'].size, grid['p1'][0], grid['p1'][-1])
            s += "      xi: {} elements from {} to {}\n".format(grid['p2'].size, grid['p2'][0], grid['p2'][-1])
        elif grid['type'] == MomentumGrid.MOMENTUMGRID_TYPE_PPARPPERP:
            s += "      ppar:  {} elements from {} to {}\n".format(grid['p1'].size, grid['p1'][0], grid['p1'][-1])
            s += "      pperp: {} elements from {} to {}\n".format(grid['p2'].size, grid['p2'][0], grid['p2'][-1])
        else:
            s += "      p1: {} elements from {} to {}\n".format(grid['p1'].size, grid['p1'][0], grid['p1'][-1])
            s += "      p2: {} elements from {} to {}\n".format(grid['p2'].size, grid['p2'][0], grid['p2'][-1])

        return s


    def setGrid(self, grid):
        """
        Set grid data based on output from DREAM.
        """
        self.t = grid['t']
        self.r = grid['r']

        # Workaround for initial data which doesn't have a time grid from DREAM
        # (TODO we should fix this in the kernel instead)
        if self.t.size == 0:
            self.t = np.array([0])

        if 'hottail' in grid:
            self.hottail = grid['hottail']

        if 'runaway' in grid:
            self.runaway = grid['runaway']

