
import matplotlib.pyplot as plt
import numpy as np
from .OtherScalarQuantity import OtherScalarQuantity


class FrozenDIp(OtherScalarQuantity):
    

    def __init__(self, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(*args, **kwargs)


    def plot(self, ax=None, show=None, **kwargs):
        """
        Override the default 'plot' routine.
        """
        genax = ax is None

        ax = super().plot(ax=ax, show=False, **kwargs)

        # Plot dIp=0
        ax.plot([min(self.time[:]), max(self.time[:])], [0, 0], 'r--')

        if genax and show is None:
            show = True

        if show:
            plt.show(block=False)
