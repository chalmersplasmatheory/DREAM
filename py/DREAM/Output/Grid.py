
import numpy as np
import DREAM.Settings.MomentumGrid as MomentumGrid
from .OutputException import OutputException
from .PXiGrid import PXiGrid
from .PparPperpGrid import PparPperpGrid


class Grid:
    

    def __init__(self, grid=None):
        """
        Constructor.

        grid: Grid data from DREAM output.
        """
        self.t = None
        self.r = None
        self.r_f = None
        self.dr = None
        self.VpVol = None
        self.hottail = None
        self.runaway = None

        # Geometric quantities
        self.effectivePassingFraction = None

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


    def getTimeAndUnit(self, t):
        """
        If 't' is an integer, returns the time at the specified index,
        together with the appropriate unit. Otherwise, if 't' is a float,
        returns the specified time with the appropriate unit.
        """
        if type(t) == int:
            return self.getTimeAndUnit(t=self.t[t])

        unit = ['s', 'ms', r'\mu s', 'ns', 'ps']
        tval = t
        ui = 0

        if tval == 0: return tval, unit[ui]

        while tval < 1 and ui+1 < len(unit):
            tval *= 1e3
            ui += 1

        return tval, unit[ui]


    def integrate(self, data, w=1.0, axis=-1):
        """
        Evaluate a numerical volume integral of the given data
        on this grid.

        data: Data to integrate.
        w:    Optional weighting function.
        axis: Axis to integrate over.
        """
        return (self.VpVol*self.dr * w * data).sum(axis)


    def strType(self, grid):
        """
        Returns the grid type as a string representation
        """
        if grid is None: return ""

        if grid.type == MomentumGrid.TYPE_PXI:
            return 'p/xi'
        elif grid.type == MomentumGrid.TYPE_PPARPPERP:
            return 'ppar/pperp'
        else:
            return '<unknown-type>'


    def strGrid(self, grid):
        """
        Returns some information about the grid as a string.
        """
        if grid is None: return ""

        s  = "      {:6s} {} elements from {} to {}\n".format(grid.p1name+':', grid.p1.size, grid.p1[0], grid.p1[-1])
        s += "      {:6s} {} elements from {} to {}\n".format(grid.p2name+':', grid.p2.size, grid.p2[0], grid.p2[-1])

        return s


    def setGrid(self, grid):
        """
        Set grid data based on output from DREAM.
        """
        self.t = grid['t']
        self.r = grid['r']
        self.r_f = grid['r_f']
        self.dr = grid['dr']
        self.VpVol = grid['VpVol']

        if 'effectivePassingFraction' in grid:
            self.effectivePassingFraction = grid['effectivePassingFraction']

        # Workaround for initial data which doesn't have a time grid from DREAM
        # (TODO we should fix this in the kernel instead)
        if self.t.size == 0:
            self.t = np.array([0])

        if 'hottail' in grid:
            self.hottail = self._generateMomentumGrid('hottail', data=grid['hottail'])

        if 'runaway' in grid:
            self.runaway = self._generateMomentumGrid('runaway', data=grid['runaway'])


    def _generateMomentumGrid(self, name, data):
        """
        Generate momentum grid object of the appropriate type.

        name: Grid name.
        data: Raw grid data from DREAM output file.
        """
        if data['type'] == MomentumGrid.TYPE_PXI:
            return PXiGrid(name=name, rgrid=self, data=data)
        elif data['type'] == MomentumGrid.TYPE_PPARPPERP:
            return PparPperpGrid(name=name, r=self, data=data)
        else:
            raise OutputException("grid: Unrecognized grid type: {}.".format(data['type']))
            

