# Special implementation for 'E_field'

import numpy as np
import scipy.constants

from . FluidQuantity import FluidQuantity
from . OutputException import OutputException


class ElectricField(FluidQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)

    
    def getNormEfield(self, field, r=None, t=None):
        """
        Returns an electric field from the other quantities by name.
        This routine is intended as a uniform interface for fetching
        quantities such as Ec, Eceff, ED etc.
        """
        # List of supported quantities (to avoid user error)
        nrm = ['Eceff', 'Ecfree', 'Ectot', 'Ec', 'ED', 'EDreic']
        if field == 'Ec': field = 'Ectot'
        elif field == 'ED': field = 'EDreic'

        if 'fluid' not in self.output.other:
            raise OutputException('No "other" fluid quantities saved in output. Normalizing electric fields are thus not available.')
        if field not in nrm:
            raise OutputException("Cannot normalize to '{}': This seems to not make sense.".format(field))
        if field not in self.output.other.fluid:
            raise OutputException("Cannot normalize to '{}': quantity not saved to output after simulation.".format(field))

        return self.output.other.fluid[field].get(r=r, t=t)


    def maxEnergy(self, t=-1):
        r"""
        Evaluates the maximum attainable runaway energy (in normalized units)
        at time ``t``. This energy is obtained by integrating the equation of
        motion:

        .. math::
        
            \frac{\mathrm{d}p}{\mathrm{d}t} = eE \quad\implies\quad
            p = \int_0^t eE(t)\,\mathrm{d}t',\\
            W = mc^2\sqrt{p^2+1},

        where :math:`e` is the elementary charge and :math:`p` is the electron
        momentum.

        :param int t: Index of time to calculate transferred momentum until.
        """
        p = self.maxMomentum(t=t)
        return np.sqrt(p**2 + 1)


    def maxMomentum(self, t=-1):
        r"""
        Evaluates the maximum attainable runaway momentum (in normalized units)
        at time ``t``. This momentum is obtained by integrating the equation of
        motion:

        .. math::
        
            \frac{\mathrm{d}p}{\mathrm{d}t} = eE \quad\implies\quad
            p = \int_0^t eE(t)\,\mathrm{d}t',

        where :math:`e` is the elementary charge and :math:`p` is the electron
        momentum.

        :param int t: Index of time to calculate transferred momentum until.
        """
        if np.isscalar(t):
            p = np.trapz(self[:t], self.grid.t[:t], axis=0)
        else:
            p = []
            t = np.asarray(t)

            if t.ndim != 1:
                raise OutputException("Unrecognized dimensions of time index: {}.".format(t.ndim))

            for time in t:
                p.append(np.trapz(self[:time], self.grid.t[:time], axis=0))

            p = np.array(p)

        p *= scipy.constants.e / (scipy.constants.m_e * scipy.constants.c)
        return p


    def norm(self, to='Ec'):
        """
        Return the value of this quantity normalized to the
        quantity specified by 'to'.

        to: Name of quantity to normalize electric field to.
            Possible values:
               Eceff   Effective critical electric field (as defined by Hesslow et al)
               Ecfree  Connor-Hastie threshold field (calculated with n=n_free)
               Ectot   Connor-Hastie threshold field (calculated with n=n_tot)
               Ec      (alias for 'Ectot')
               ED      Dreicer field
            Note that the quantity with which to normalize to
            must be saved as an 'OtherQuantity'.
        """
        return self.normalize(to=to)


    def normalize(self, to='Ec'):
        """
        Return the value of this quantity normalized to the
        quantity specified by 'to'.

        to: Name of quantity to normalize electric field to.
            Possible values:
               Eceff   Effective critical electric field (as defined by Hesslow et al)
               Ecfree  Connor-Hastie threshold field (calculated with n=n_free)
               Ectot   Connor-Hastie threshold field (calculated with n=n_tot)
               Ec      (alias for 'Ectot')
               EDreic  Dreicer field
               ED      (alias for 'EDreic')
            Note that the quantity with which to normalize to
            must be saved as an 'OtherQuantity'.
        """
        # OtherQuantity's are not defined at t=0, so we extend them
        # arbitrarily here (in order for the resulting FluidQuantity to
        # be plottable on the regular time grid)
        Enorm = np.zeros(self.data.shape)
        Enorm[1:,:] = self.getNormEfield(field=to)
        Enorm[0,:]  = Enorm[1,:]

        data = self.data[:] / Enorm
        return FluidQuantity(name='E / {}'.format(to), data=data, grid=self.grid, output=self.output)


    def plot(self, norm=None, **kwargs):
        """
        Wrapper for FluidQuantity.plot(), adding the 'norm' argument
        (for directly plotting a normalized electric field)
        """
        if norm is None:
            return super().plot(**kwargs)
        else:
            _E = self.normalize(to=norm)
            return _E.plot(**kwargs)


    def plotRadialProfile(self, norm=None, **kwargs):
        """
        Wrapper for FluidQuantity.plotRadialProfile(), adding the 'norm'
        argument (for directly plotting a normalized electric field)
        """
        if norm is None:
            return super().plotRadialProfile(**kwargs)
        else:
            _E = self.normalize(to=norm)
            return _E.plotRadialProfile(**kwargs)


    def plotTimeProfile(self, norm=None, **kwargs):
        """
        Wrapper for FluidQuantity.plotTimeProfile(), adding the 'norm'
        argument (for directly plotting a normalized electric field)
        """
        if norm is None:
            return super().plotTimeProfile(**kwargs)
        else:
            _E = self.normalize(to=norm)
            return _E.plotTimeProfile(**kwargs)


    def plotIntegral(self, norm=None, **kwargs):
        """
        Wrapper for FluidQuantity.plotTimeProfile(), adding the 'norm'
        argument (for directly plotting a normalized electric field)
        """
        if norm is None:
            return super().plotIntegral(**kwargs)
        else:
            _E = self.normalize(to=norm)
            return _E.plotIntegral(**kwargs)


