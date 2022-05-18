

import matplotlib.pyplot as plt
import numpy as np
from . ScalarQuantity import ScalarQuantity


class PlasmaCurrent(ScalarQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.

        :param str name:           Name of quantity.
        :param numpy.ndarray data: Data for quantity.
        :param list attr:          List of attributes for quantity.
        :param grid:               Grid on which the quantity is defined.
        :param DREAMOutput output: Parent output object.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)


    def plotCurrents(self, ax=None, show=None, t=None):
        """
        Plot the runaway, hot, and ohmic currents, alongside the plasma
        current in the same plot.
        """
        Ire  = self.output.eqsys.j_re.current()
        Ihot = self.output.eqsys.j_hot.current()
        Iohm = self.output.eqsys.j_ohm.current()

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        if t is None:
            t = slice(None)

        ax.plot(self.time[t], self.data[t], 'k', label=r'$I_{\rm p}$')
        ax.plot(self.time[t], Ire[t], 'r--', label=r'$I_{\rm re}$')
        ax.plot(self.time[t], Iohm[t], 'b--', label=r'$I_{\Omega}$')
        ax.plot(self.time[t], Ihot[t], 'y--', label=r'$I_{\rm hot}$')

        ax.set_xlabel(r'Time $t$ (s)')
        ax.set_ylabel(r'Plasma current')

        ax.set_xlim([min(self.time[t]), max(self.time[t])])

        Imax = max(np.amax(abs(self.data[t])), np.amax(abs(Ire)), np.amax(abs(Ihot)), np.amax(abs(Iohm)))
        if np.any(self.data[t]<0) or np.any(Ire<0) or np.any(Ihot<0) or np.any(Iohm<0):
            ax.set_ylim([-1.1*Imax, 1.1*Imax])
        else:
            ax.set_ylim([0, 1.1*Imax])
        
        ax.legend(frameon=False)

        if show:
            plt.show(block=False)

        return ax


