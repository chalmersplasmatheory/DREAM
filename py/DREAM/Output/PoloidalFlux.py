# Special treatment for poloidal flux

import numpy as np

from .FluidQuantity import FluidQuantity
from .OutputException import OutputException


class PoloidalFlux(FluidQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)


    def getNormalized(self, t=None):
        """
        Returns the normalized poloidal flux at the specified time.
        This poloidal flux could be used as a radial coordinate.
        """
        if t is None:
            t = slice(None)

        # Extrapolate psi_p to get psi_p(r=0)
        dr = self.grid.r[0] / (self.grid.r[1]-self.grid.r[0])
        psi0 = self.data[t,0]*(1+dr) - self.data[t,1]*dr
        psi_edge = self.output.eqsys.psi_edge[t,0]

        psi_n = (self.data[t,:] - psi0) / (psi_edge - psi0)

        return psi_n


