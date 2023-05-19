
import matplotlib.pyplot as plt
import numpy as np
from .FluidQuantity import FluidQuantity
from .OutputException import OutputException


class RunawayElectronDensity(FluidQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, grid=grid, output=output, attr=attr)


    def plotRates(self, r=None, t=None, ax=None, show=True):
        """
        Plot all runaway rates in the same figure.
        """
        if t is None and r is None:
            #raise OutputException("When plotting all runaway rates, at least one of 'r' and 't' must be specified.")
            integrate = True

        labels = []

        if 'fluid' not in self.output.other:
            raise OutputException("No 'other' fluid quantities saved in output. Cannot plot runaway rates.")

        # Plot total runaway rate
        if 'runawayRate' in self.output.other.fluid:
            if integrate:
                ax = self.output.other.fluid.runawayRate.plotIntegral(ax=ax, show=False)
            else:
                ax = self.output.other.fluid.runawayRate.plot(r=r, t=t, ax=ax, show=False)
            labels.append('Total')

        # Plot avalanche
        if 'GammaAva' in self.output.other.fluid:
            ax = self.output.other.fluid.GammaAva.plotRunawayRate(r=r, t=t, ax=ax, show=False)
            labels.append('Avalanche')

        for o in self.output.other.fluid.keys():
            if not o.startswith('gamma'):
                continue

            q = self.output.other.fluid[o]
            
            # Ignore if empty
            if np.sum(np.abs(q[:])) == 0:
                continue

            if integrate:
                ax = q.plotIntegral(ax=ax, show=False)
            else:
                ax = q.plot(r=r, t=t, ax=ax, show=False)
            labels.append(o[5:].replace(r'_', r'\_'))

        plt.ylabel('Runaway rates')
        plt.legend(labels)

        if show:
            plt.show(block=False)

        return ax

