# ppar/pperp output momentum grid object


import numpy as np
from .MomentumGrid import MomentumGrid


class PparPperpGrid(MomentumGrid):
    

    def __init__(self, name, r, dr, data):
        """
        Constructor.

        name: Grid name.
        r:    Associated radial grid.
        data: Momentum grid data.
        """
        super(PXiGrid, self).__init__(name=name, r=r, dr=dr, data=data)

        self.ppar   = data['p1']
        self.pperp  = data['p2']
        self.dppar  = data['dp1']
        self.dpperp = data['dp2']

