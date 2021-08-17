# Special quantity for ion thermal energy


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

