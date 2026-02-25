# Special quantity for ion thermal energy


import matplotlib.pyplot as plt
import scipy.constants
from . FluidQuantity import FluidQuantity
from . IonSpeciesFluidQuantity import IonSpeciesFluidQuantity


class IonThermalEnergy(IonSpeciesFluidQuantity):
    

    def __init__(self, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(*args, **kwargs)


    def getTemperature(self, ion=None):
        """
        Returns the temperature of the named ion species. If no
        ion name is provided, a compound ion temperature object
        is returned, containing T_i for all ions.
        """
        ec = scipy.constants.e

        if ion is None:
            W_i = self.data[:]
            N_i = self.output.eqsys.N_i.data[:]
            T_i = (2/3) * (W_i/ec) / N_i

            return IonSpeciesFluidQuantity(name='T_i', data=T_i, grid=self.grid, output=self.output, attr=self.attr)
        else:
            W_i = self[ion][:]
            N_i = self.output.eqsys.N_i[ion][:]
            T_i = (2/3) * (W_i/ec) / N_i

            return FluidQuantity(name=f'Ti_{ion}', data=T_i, grid=self.grid, output=self.output, attr=self.attr)


    def plotTemperature(self, ion=None, ax=None, show=None, r=None, t=None, *args, **kwargs):
        """
        Plot ion temperature
        """
        T_i = self.getTemperature(ion)
        return T_i.plot(ax=ax, show=show, r=r, t=t, *args, **kwargs)


    def plotEnergyBalance(self, r=None, t=None, ax=None, show=True, log=False):
        """
        Plot all the terms appearing in the ion energy balance equation.
        Plots quantities for each ion species.
        """
        integrate = False
        if t is None and r is None:
            integrate = True

        labels = []
        for o in self.output.other.fluid.keys():
            if not o.startswith('Ti_'):
                continue

            q = self.output.other.fluid[o]
            
            for ion_name in q.ions.getNames():
                ion_q = q[ion_name]
                
                if integrate:
                    ax = ion_q.plotIntegral(ax=ax, show=False)
                else:
                    ax = ion_q.plot(r=r, t=t, ax=ax, show=False, log=log)
                
                labels.append('{} ({})'.format(o[3:].replace(r'_', r'\_'), ion_name))

        plt.legend(labels)

        if show:
            plt.ylabel(r'$\dot{W}_i\ \mathrm{(W)}$')
            plt.show(block=False)

        return ax

