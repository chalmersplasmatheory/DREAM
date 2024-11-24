
import matplotlib.pyplot as plt
from .OtherFluidQuantity import OtherFluidQuantity
from .OutputException import OutputException


class AvalancheGrowthRate(OtherFluidQuantity):
    

    def __init__(self, name, data, description, grid, output, momentumgrid):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, description=description, grid=grid, output=output)


    def plotNormalized(self, r=None, t=None, ax=None, show=True, norm='Eceff'):
        """
        Plot the avalanche growth rate, normalized to E-Eceff.
        """
        ut = self._renormalizeTimeIndexForUnknown(t)

        Enorm = self.output.eqsys.E_field.getNormEfield(field=norm, r=r, t=t)
        E     = self.output.eqsys.E_field.get(r=r, t=ut)
        EE    = E-Enorm

        return self.plot(r=r, t=t, ax=ax, show=show, weight=EE)


    def plotRunawayRate(self, r=None, t=None, ax=None, show=True):
        """
        Plot the runaway corresponding to this growth rate (i.e.
        GammaAva*n_re).
        """
        ut = self._renormalizeTimeIndexForUnknown(t)

        n_re = self.output.eqsys.n_re.get(r=r, t=ut)
        if t is None and r is None:
            return self.plotIntegral(ax=ax, show=show, w=n_re)
        else:
            return self.plot(r=r, t=t, ax=ax, show=show, weight=n_re)

              
    def getRunawayRate(self, r=None, t=None, ax=None, show=True):
        """
        Calculates the runaway rate corresponding to this growth rate (i.e.
        GammaAva*n_re).
        """
        ut = self._renormalizeTimeIndexForUnknown(t)

        n_re = self.output.eqsys.n_re.get(r=r, t=ut)
        return self.data*n_re
        
        
    def getRunawayRateIntegral(self, r=None, t=None, ax=None, show=True):
        """
        Calculates the runaway rate corresponding to this growth rate (i.e.
        GammaAva*n_re).
        """
        ut = self._renormalizeTimeIndexForUnknown(t)

        n_re = self.output.eqsys.n_re.get(r=r, t=ut)
        return self.integral(w=n_re)
        
