# Class for managing diagnostics about the equation trigger

import matplotlib.pyplot as plt
import numpy as np


class TriggerInformation:
    

    def __init__(self, qty, data, output):
        """
        Constructor.
        """
        self.name = qty.name
        self.grid = qty.grid
        self.data = data
        self.output = output


    def __getitem__(self, idx):
        return self.data[idx]


    def plot(self, ax=None, show=None, r=None, t=None, **kwargs):
        """
        Visualize the trigger state.

        :param ax:       Matplotlib axes object to use for plotting.
        :param show:     If ``True``, shows the plot immediately via a call to ``matplotlib.pyplot.show()`` with ``block=False``. If ``None``, this is interpreted as ``True`` if ``ax`` is also ``None``.
        :param r:        Radial index/ices to show. If ``None``, shows all radial points.
        :param t:        Time index/ices to show. If ``None``, shows all time points.
        :param **kwargs: Arguments passed to ``matplotlib.Axes.pcolor()``.

        :return: a matplotlib axis object.
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        if r is None:
            r = slice(None)
        if t is None:
            t = slice(None)

        data = self.data[t,r]
        nr = self.data.shape[-1]
        rarr = np.linspace(0, nr-1, nr)[r]
        tarr = self.grid.t[t]

        ax.pcolormesh(rarr, tarr, data, cmap='RdYlGn', shading='nearest', vmin=0, vmax=1)
        ax.set_xlabel('Index')
        ax.set_ylabel('Time (s)')

        if show:
            plt.show(block=False)

        return ax


