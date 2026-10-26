"""Neutral thermal energy density [J/m^3] and temperature [eV]."""

import numpy as np
import scipy.constants
from .FluidQuantity import FluidQuantity
from .IonSpeciesFluidQuantity import IonSpeciesFluidQuantity


class NeutralThermalEnergy(IonSpeciesFluidQuantity):
    def getTemperature(self, ion=None):
        """Return neutral temperatures, using only charge-state-zero densities.

        Temperature is NaN where the neutral density is zero or negative.
        The factor 3/2 matches the W_n definition in the DREAM solver.
        """
        if ion is None:
            temperatures = []
            for name in self.ions.names:
                temperatures.append(self.getTemperature(name).data)

            return IonSpeciesFluidQuantity(
                name='T_n', data=np.stack(temperatures, axis=1),
                grid=self.grid, output=self.output, attr=self.attr
            )

        energy = self[ion][:]
        density = self.output.eqsys.n_i[ion][0][:]
        temperature = np.full(energy.shape, np.nan, dtype=float)
        np.divide(
            energy, 1.5 * scipy.constants.e * density,
            out=temperature, where=density > 0
        )

        return FluidQuantity(
            name=f'Tn_{ion}', data=temperature,
            grid=self.grid, output=self.output, attr=self.attr
        )

    def plotTemperature(self, ion=None, ax=None, show=None, r=None, t=None,
                        *args, **kwargs):
        """Plot neutral temperatures with the same interface as W_i."""
        temperature = self.getTemperature(ion)
        return temperature.plot(
            ax=ax, show=show, r=r, t=t, *args, **kwargs
        )
