# Settings object for the hot-tail grid.

import numpy as np


from PGrid import PGrid
from XiGrid import XiGrid


class HotTailGrid():

    MOMENTUMGRID_TYPE_PXI = 1
    MOMENTUMGRID_TYPE_PPARPPERP = 2
    
    def __init__(self, enabled=True, ttype=1, np=100, nxi=1, pmax=None, pgrid=None, xigrid=None):
        """
        Constructor.

        enabled: If 'True', enables the hot-tail grid in the simulation.
        ttype:   Type of momentum grid (p/xi or ppar/pperp).
        np:      Number of momentum grid points.
        nxi:     Number of pitch grid points.
        pmax:    Maximum momentum on grid.
        pgrid:   Momentum grid object.
        xigrid:  Pitch grid object.
        """
        if pgrid is None: pgrid = PGrid()
        if xigrid is None: xigrid = XiGrid()

        self.set(enabled=enabled, ttype=ttype, np=np, nxi=nxi, pmax=pmax, pgrid=pgrid, xigrid=xigrid)


    def set(self, enabled=True, ttype=1, np=100, nxi=1, pmax=None):
        """
        Set all settings for this hot-tail grid.
        """
        self.enabled = enabled

        self.ttype = ttype
        if self.ttype == self.MOMENTUMGRID_TYPE_PXI:
            self.pgrid = PGrid()
            self.xigrid = XiGrid()

            self.setNp(np)
            self.setNxi(nxi)
            self.setPmax(pmax)
        elif self.ttype == self.MOMENTUMGRID_TYPE_PPARPPERP:
            raise DREAMException("No support implemented yet for 'ppar/pperp' grids.")
        else:
            raise DREAMException("Unrecognized momentum grid type specified: {}.".format(ttype))

    ##################
    # SETTERS
    ##################
    def setNp(self, np):
        if np <= 0:
            raise DREAMException("HotTailGrid: Invalid value assigned to 'np': {}. Must be > 0.".format(np))
        elif np == 1:
            print("WARNING: HotTailGrid: np = 1. Consider disabling the hot-tail grid altogether.")

        self.np = int(np)
        self.pgrid.setNp(self.np)


    def setNxi(self, nxi):
        if nxi <= 0:
            raise DREAMException("HotTailGrid: Invalid value assigned to 'nxi': {}. Must be > 0.".format(nxi))

        self.nxi = int(nxi)
        self.xigrid.setNp(self.np)


    def setPmax(self, pmax):
        if pmax <= 0:
            raise DREAMException("HotTailGrid: Invalid value assigned to 'pmax': {}. Must be > 0.".format(pmax))

        self.pmax = float(pmax)


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this HotTailGrid object.
        """
        data = {
            'enabled': self.enabled,
            'type':    self.type
        }

        if self.type == MOMENTUMGRID_TYPE_PXI:
            data = {**data, **(self.pgrid.todict()), **(self.xigrid.todict())}
        elif self.type == MOMENTUMGRID_TYPE_PPARPPERP:
            raise DREAMException("No support implemented yet for 'ppar/pperp' grids.")
        else:
            raise DREAMException("Unrecognized momentum grid type specified: {}.".format(ttype))

        return data


