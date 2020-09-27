# Object representing a single ion charge state


import numpy as np
from .FluidQuantity import FluidQuantity


class IonState(FluidQuantity):
    

    def __init__(self, name, Z, Z0, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=("{}-{}".format(name, Z0)), data=data, attr=attr, grid=grid, output=output)

        self.Z  = Z
        self.Z0 = Z0

