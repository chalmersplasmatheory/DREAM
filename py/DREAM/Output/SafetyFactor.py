
import matplotlib.pyplot as plt
from .OtherFluidQuantity import OtherFluidQuantity
from .OutputException import OutputException


class SafetyFactor(OtherFluidQuantity):
    

    def __init__(self, name, data, description, grid, output, momentumgrid):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, description=description, grid=grid, output=output)


    def plotq(self, r=None, t=None, ax=None, show=True):
        """
        Plot only the safety factor (instead of ``q*R0``).
        """
        return self.plot(r=r, t=t, ax=ax, show=show, weight=1/self.grid.R0[0])


