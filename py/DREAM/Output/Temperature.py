

import matplotlib.pyplot as plt
from .FluidQuantity import FluidQuantity
from .OutputException import OutputException


class Temperature(FluidQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, grid=grid, output=output, attr=attr)


    def plotEnergyBalance(self, r=None, t=None, ax=None, show=True, log=False):
        """
        Plot all the terms appearing in the energy balance equation.
        """
        integrate = False
        if t is None and r is None:
            integrate = True

        labels = []
        for o in self.output.other.fluid.keys():
            if not o.startswith('Tcold_'):
                continue

            q = self.output.other.fluid[o]

            if integrate:
                ax = q.plotIntegral(ax=ax, show=show)
            else:
                ax = q.plot(r=r, t=t, ax=ax, show=False, log=log)

            labels.append(o[6:].replace(r'_', r'\_'))

        if integrate and 'energyloss_T_cold' in self.output.other.scalar.keys():
            o = self.output.other.scalar.energyloss_T_cold
            ax = o.plot(ax=ax, show=show)
            labels.append(o.name.replace(r'_', r'\_'))

        plt.legend(labels)

        if show:
            plt.show(block=False)

        return ax


