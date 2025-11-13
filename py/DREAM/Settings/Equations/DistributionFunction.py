
import numpy as np
from . EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
from .. import AdvectionInterpolation
from .. TransportSettings import TransportSettings


# BOUNDARY CONDITIONS (WHEN f_re IS DISABLED)
BC_F_0        = 1
BC_PHI_CONST  = 2
BC_DPHI_CONST = 3

# Interpolation methods for advection term in kinetic equation
# (we keep these for backwards compatibility)
AD_INTERP_CENTRED  = AdvectionInterpolation.AD_INTERP_CENTRED
AD_INTERP_UPWIND   = AdvectionInterpolation.AD_INTERP_UPWIND
AD_INTERP_UPWIND_2ND_ORDER = AdvectionInterpolation.AD_INTERP_UPWIND_2ND_ORDER
AD_INTERP_DOWNWIND = AdvectionInterpolation.AD_INTERP_DOWNWIND
AD_INTERP_QUICK    = AdvectionInterpolation.AD_INTERP_QUICK 
AD_INTERP_SMART    = AdvectionInterpolation.AD_INTERP_SMART 
AD_INTERP_MUSCL    = AdvectionInterpolation.AD_INTERP_MUSCL 
AD_INTERP_OSPRE    = AdvectionInterpolation.AD_INTERP_OSPRE 
AD_INTERP_TCDF     = AdvectionInterpolation.AD_INTERP_TCDF  

AD_INTERP_JACOBIAN_LINEAR = AdvectionInterpolation.AD_INTERP_JACOBIAN_LINEAR
AD_INTERP_JACOBIAN_FULL   = AdvectionInterpolation.AD_INTERP_JACOBIAN_FULL  
AD_INTERP_JACOBIAN_UPWIND = AdvectionInterpolation.AD_INTERP_JACOBIAN_UPWIND

SYNCHROTRON_MODE_NEGLECT = 1
SYNCHROTRON_MODE_INCLUDE = 2

RIPPLE_MODE_NEGLECT = 1
RIPPLE_MODE_BOX = 2
RIPPLE_MODE_GAUSSIAN = 3

TIME_VARYING_B_MODE_NEGLECT = 1
TIME_VARYING_B_MODE_INCLUDE = 2

DISTRIBUTION_MODE_NUMERICAL = 1
DISTRIBUTION_MODE_ANALYTICAL = 2
DISTRIBUTION_MODE_PRESCRIBED = 3

class DistributionFunction(UnknownQuantity):
    

    def __init__(self, settings, name, grid,
        f=[0], initr=[0], initp=[0], initxi=[0],
        initppar=None, initpperp=None,
        rn0=None, n0=None, rT0=None, T0=None, bc=BC_PHI_CONST,
        ad_int_r=AD_INTERP_CENTRED, ad_int_p1=AD_INTERP_CENTRED,
        ad_int_p2=AD_INTERP_CENTRED, ad_jac_r=AD_INTERP_JACOBIAN_FULL,
        ad_jac_p1=AD_INTERP_JACOBIAN_FULL, ad_jac_p2=AD_INTERP_JACOBIAN_FULL,
        mode = DISTRIBUTION_MODE_NUMERICAL, fluxlimiterdamping=1.0):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.name = name
        self.grid = grid

        self.boundarycondition = bc
        
        self.mode = mode
        self.ripplemode = RIPPLE_MODE_NEGLECT
        self.synchrotronmode = SYNCHROTRON_MODE_NEGLECT
        self.timevaryingbmode = TIME_VARYING_B_MODE_NEGLECT
        self.transport = TransportSettings(kinetic=True)
        self.fullIonJacobian = True

        self.advectionInterpolation = AdvectionInterpolation.AdvectionInterpolation(
            kinetic=True,
            ad_int_r=ad_int_r, ad_int_p1=ad_int_p1, ad_int_p2=ad_int_p2,
            ad_jac_r=ad_jac_r, ad_jac_p1=ad_jac_p1, ad_jac_p2=ad_jac_p2,
            fluxlimiterdamping=fluxlimiterdamping)

        self.n0  = rn0
        self.rn0 = n0

        self.T0  = rT0
        self.rT0 = T0

        self.prescribed_t = None
        self.prescribed_r = None
        self.prescribed_f = None
        self.prescribed_p = None
        self.prescribed_xi = None
        self.prescribed_ppar = None
        self.prescribed_pperp = None

        self.init = None

        if f is not None:
            self.setInitialValue(f, r=initr, p=initp, xi=initxi, ppar=initppar, pperp=initpperp)
        elif n0 is not None:
            self.setInitialProfiles(rn0=rn0, n0=n0, rT0=rT0, T0=T0)


    def setBoundaryCondition(self, bc):
        """
        Sets the boundary condition at p=pmax. For 'f_hot', this boundary
        condition is only used when 'f_re' is disabled.

        :param int bc: Flag specifying which boundary condition to use.
        """
        self.boundarycondition = bc

    def setAdvectionInterpolationMethod(self,ad_int=None, ad_int_r=AD_INTERP_CENTRED,
        ad_int_p1=AD_INTERP_CENTRED, ad_int_p2=AD_INTERP_CENTRED, ad_jac=None, 
        ad_jac_r=AD_INTERP_JACOBIAN_FULL, ad_jac_p1=AD_INTERP_JACOBIAN_FULL,
        ad_jac_p2=AD_INTERP_JACOBIAN_FULL, fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used in the advection terms of
        the kinetic equation. To set all three components, provide ad_int and/or ad_jac.
        Otherwise the three components can use separate interpolation methods.
        
        :param int ad_int:               Interpolation method to use for all coordinates.
        :param int ad_int_r:             Interpolation method to use for the radial coordinate.
        :param int ad_int_p1:            Interpolation method to use for the first momentum coordinate.
        :param int ad_int_p2:            Interpolation method to use for the second momentum coordinate.
        :param int ad_jac:               Jacobian interpolation mode to use for all coordinates.
        :param int ad_jac_r:             Jacobian interpolation mode to use for the radial coordinate.
        :param int ad_jac_p1:            Jacobian interpolation mode to use for the first momentum coordinate.
        :param int ad_jac_p2:            Jacobian interpolation mode to use for the second momentum coordinate.
        :param float fluxlimiterdamping: Damping parameter used to under-relax the interpolation coefficients during non-linear iterations (should be between 0 and 1).
        """
        self.advectionInterpolation.setMethod(ad_int=ad_int, ad_int_r=ad_int_r,
            ad_int_p1=ad_int_p1, ad_int_p2=ad_int_p2, ad_jac=ad_jac, 
            ad_jac_r=ad_jac_r, ad_jac_p1=ad_jac_p1,
            ad_jac_p2=ad_jac_p2, fluxlimiterdamping=fluxlimiterdamping)


    def setInitialProfiles(self, n0, T0, rn0=None, rT0=None):
        """
        Sets the initial density and temperature profiles of the electron
        population.

        :param rn0: Radial grid on which the density is given.
        :param n0:  Electron density profile.
        :param rT0: Radial grid on which the temperature is given.
        :param T0:  Electron temperature profile.
        """
        if rn0 is not None:
            self.rn0 = np.asarray(rn0)
        else:
            if not np.isscalar(n0):
                raise EquationException("{}: Non-scalar initial density profile given, but no radial grid specified.".format(self.name))
            self.rn0 = np.array([0])

        if rT0 is not None:
            self.rT0 = np.asarray(rT0)
        else:
            if not np.isscalar(T0):
                raise EquationException("{}: Non-scalar initial temperature profile given, but no radial grid specified.".format(self.name))
            self.rT0 = np.array([0])

        self.n0  = np.asarray(n0)
        self.T0  = np.asarray(T0)

        if self.rn0.ndim == 0: self.rn0 = np.asarray([self.rn0])
        if self.n0.ndim == 0:  self.n0 = np.asarray([self.n0])
        if self.rT0.ndim == 0: self.rT0 = np.asarray([self.rT0])
        if self.T0.ndim == 0:  self.T0 = np.asarray([self.T0])

        # Reset numerically provided distribution (if any)
        self.init = None

        self.verifyInitialProfiles()


    def setInitialValue(self, f, r, p=None, xi=None, ppar=None, pperp=None):
        """
        Set the initial value of this electron distribution function. Only one
        of the pairs (p, xi) and (ppar, pperp) of momentum grids need to be
        given.

        :param f:     Array representing the distribution function value on the grid (must have size (nr, nxi, np) or (nr, npperp, nppar))
        :param r:     Radial grid on which the initial distribution is given.
        :param p:     Momentum grid.
        :param xi:    Pitch grid.
        :param ppar:  Parallel momentum grid.
        :param pperp: Perpendicular momentum grid.
        """
        self.init = {}

        def conv(v):
            if type(v) == list:
                return np.array(v)
            elif type(v) == float or type(v) == int:
                return np.array([float(v)])
            else:
                return v

        ff = conv(f)
        self.init['r'] = conv(r)

        if p is not None and xi is not None:
            self.init['p'] = conv(p)
            self.init['xi'] = conv(xi)
            self.init['ppar'] = np.array([])
            self.init['pperp'] = np.array([])

            if ff.size == 1:
                ff = ff * np.ones((self.init['r'].size, self.init['xi'].size, self.init['p'].size))
        elif ppar is not None and pperp is not None:
            self.init['ppar'] = conv(ppar)
            self.init['pperp'] = conv(pperp)
            self.init['p'] = np.array([])
            self.init['xi'] = np.array([])

            if ff.size == 1:
                ff = ff * np.ones((self.init['r'].size, self.init['pperp'].size, self.init['ppar'].size))
        else:
            raise EquationException("{}: No momentum grid given for initial value.".format(self.name))

        self.init['x'] = ff

        # Reset initial profiles (if any)
        self.rn0 = self.rT0 = None
        self.n0 = self.T0 = None

        self.verifyInitialDistribution()


    def enableAnalyticalDistribution(self, mode=True):
        """
        Enables/disables the use of an analytical distribution
        function to represent the electron population
        """
        if mode:
            self.mode = DISTRIBUTION_MODE_ANALYTICAL
        else:
            self.mode = DISTRIBUTION_MODE_NUMERICAL


    def setRippleMode(self, mode):
        """
        Enables/disables inclusion of pitch scattering due to the magnetic ripple.

        :param int mode: Flag indicating whether or not to include magnetic ripple effects.
        """
        if type(mode) == bool:
            self.ripplemode = RIPPLE_MODE_BOX if mode else RIPPLE_MODE_NEGLECT
        else:
            self.ripplemode = int(mode)


    def setTimeVaryingB(self, mode):
        """
        Enables/disable the time-varying magnetic field strength operator.

        :param int mode: Flag indicating whether or not to include the time-varying magnetic field operator.
        """
        if type(mode) == bool:
            self.timevaryingbmode = TIME_VARYING_B_MODE_INCLUDE if mode else TIME_VARYING_B_MODE_NEGLECT
        else:
            self.timevaryingbmode = int(mode)


    def setSynchrotronMode(self, mode):
        """
        Sets the type of synchrotron losses to have (either enabled or disabled).

        :param int mode: Flag indicating whether or not to enable synchrotron losses (may be bool).
        """
        if type(mode) == bool:
            self.synchrotronmode = SYNCHROTRON_MODE_INCLUDE if mode else SYNCHROTRON_MODE_NEGLECT
        else:
            self.synchrotronmode = int(mode)

    def enableIonJacobian(self, includeJacobian):
        """
        Enables/disables the ion jacobian in the kinetic equation.

        :param bool includeJacobian: Flag indicating whether the ion jacobian will be added. True by default, False to disable.
        """
        self.fullIonJacobian = includeJacobian


    def prescribe(self, f, t, r, xi=None, p=None, pperp=None, ppar=None):
        """
        Prescribe the time evolution of this distribution function instead of
        solving for it using a kinetic equation.
        """
        f = np.array(f)
        t = np.array(t)
        r = np.array(r)
        if p is not None and xi is not None:
            p = np.array(p)
            xi = np.array(xi)
        elif ppar is not None and pperp is not None:
            ppar = np.array(ppar)
            pperp = np.array(pperp)
        else:
            raise EquationException("Either 'p' and 'xi', or 'ppar' and 'pperp', must be specified.")

        if t.ndim != 1:
            raise EquationException("The time vector 't' must be a one-dimensional array.")
        if r.ndim != 1:
            raise EquationException("The radius vector 'r' must be a one-dimensional array.")

        self.mode = DISTRIBUTION_MODE_PRESCRIBED 
        self.prescribed_t = t
        self.prescribed_r = r
        self.prescribed_f = f

        if p is not None:
            if xi.ndim != 1:
                raise EquationException("The pitch vector 'xi' must be a one-dimensional array.")
            if p.ndim != 1:
                raise EquationException("The momentum vector 'p' must be a one-dimensional array.")

            if f.shape != (t.size, r.size, xi.size, p.size):
                raise EquationException(f"The distribution function 'f' must have shape (nt, nr, nxi, np) = ({t.size}, {r.size}, {xi.size}, {p.size})")

            self.prescribed_p = p
            self.prescribed_xi = xi
        else:
            if pperp.ndim != 1:
                raise EquationException("The perpendicular momentum vector 'pperp' must be a one-dimensional array.")
            if ppar.ndim != 1:
                raise EquationException("The parallel momentum vector 'ppar' must be a one-dimensional array.")

            if f.shape != (t.size, r.size, ppar.size, pperp.size):
                raise EquationException(f"The distribution function 'f' must have shape (nt, nr, npperp, nppar) = ({t.size}, {r.size}, {pperp.size}, {ppar.size})")

            self.prescribed_ppar = None
            self.prescribed_pperp = None


    def fromdict(self, data):
        """
        Load data for this object from the given dictionary.

        :param dict data: Dictionary to load distribution function from.
        """
        def scal(v):
            if type(v) == np.ndarray: return v[0]
            else: return v

        if 'mode' in data:
            self.mode = data['mode']
        if 'boundarycondition' in data:
            self.boundarycondition = data['boundarycondition']

        if 'adv_interp' in data:
            self.advectionInterpolation.fromdict(data['adv_interp'])

        if 'init' in data:
            self.init = data['init']
        elif ('n0' in data) and ('T0' in data):
            self.rn0 = data['n0']['r']
            self.n0  = data['n0']['x']
            self.rT0 = data['T0']['r']
            self.T0  = data['T0']['x']

        if 'ripplemode' in data:
            self.ripplemode = int(scal(data['ripplemode']))

        if 'timevaryingbmode' in data:
            self.timevaryingbmode = int(scal(data['timevaryingbmode']))

        if 'synchrotronmode' in data:
            self.synchrotronmode = data['synchrotronmode']
            if type(self.synchrotronmode) != int:
                self.synchrotronmode = int(self.synchrotronmode[0])

        if 'transport' in data:
            self.transport.fromdict(data['transport'])

        if 'fullIonJacobian' in data:
            self.fullIonJacobian = bool(data['fullIonJacobian'])

        if 'f_prescribed' in data:
            self.prescribed_t = data['f_prescribed']['t']
            self.prescribed_r = data['f_prescribed']['r']
            self.prescribed_f = data['f_prescribed']['x']

            if 'p' in data['f_prescribed']:
                self.prescribed_p = data['f_prescribed']['p']
                self.prescribed_xi = data['f_prescribed']['xi']
            elif 'ppar' in data['f_prescribed']:
                self.prescribed_ppar = data['f_prescribed']['ppar']
                self.prescribed_pperp = data['f_prescribed']['pperp']
            else:
                raise EquationException("Expected either 'p' and 'xi', or 'ppar' and 'pperp', to be present under 'f_prescribed'.")

        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of this
        DistributionFunction object.

        :return: a dictionary, containing all settings of this object, which can be directly given to DREAM.
        """
        data = {}
        data['mode'] = self.mode
        if self.grid.enabled:
            data['boundarycondition'] = self.boundarycondition

            # Advection interpolation
            data['adv_interp'] = self.advectionInterpolation.todict()

            if self.init:
                data['init'] = {}
                data['init']['x'] = self.init['x']
                data['init']['r'] = self.init['r']

                if self.init['p'].size > 0 and self.init['xi'].size > 0:
                    data['init']['p'] = self.init['p']
                    data['init']['xi'] = self.init['xi']
                elif self.init['ppar'].size > 0 and self.init['pperp'].size > 0:
                    data['init']['ppar'] = self.init['ppar']
                    data['init']['pperp'] = self.init['pperp']
            elif self.n0 is not None:
                data['n0'] = { 'r': self.rn0, 'x': self.n0 }
                data['T0'] = { 'r': self.rT0, 'x': self.T0 }
            
            data['ripplemode'] = self.ripplemode
            data['synchrotronmode'] = self.synchrotronmode
            data['timevaryingbmode'] = self.timevaryingbmode
            data['transport'] = self.transport.todict()
            data['fullIonJacobian'] = self.fullIonJacobian

            if self.mode == DISTRIBUTION_MODE_PRESCRIBED:
                data['f_prescribed'] = {
                    't': self.prescribed_t,
                    'r': self.prescribed_r,
                    'x': self.prescribed_f
                }

                if self.prescribed_p is not None:
                    data['f_prescribed']['p'] = self.prescribed_p
                    data['f_prescribed']['xi'] = self.prescribed_xi
                else:
                    data['f_prescribed']['ppar'] = self.prescribed_ppar
                    data['f_prescribed']['pperp'] = self.prescribed_pperp

        if self.mode == DISTRIBUTION_MODE_ANALYTICAL:
            data['n0'] = { 'r': self.rn0, 'x': self.n0 }
            data['T0'] = { 'r': self.rT0, 'x': self.T0 }

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.grid.enabled:
            if self.mode not in [DISTRIBUTION_MODE_NUMERICAL, DISTRIBUTION_MODE_PRESCRIBED]:
                raise EquationException("{}: Invalid mode set. Must be 'NUMERICAL' or 'PRESCRIBED' when the grid is 'enabled'.".format(self.name))
            bc = self.boundarycondition
            if (bc != BC_F_0) and (bc != BC_PHI_CONST) and (bc != BC_DPHI_CONST):
                raise EquationException("{}: Invalid external boundary condition set: {}.".format(self.name, bc))
            if self.init is not None:
                self.verifyInitialDistribution()
            elif (self.n0 is not None) or (self.T0 is not None):
                self.verifyInitialProfiles()
            else:
                raise EquationException("{}: Invalid/no initial condition set for the distribution function.".format(self.name))

            self.advectionInterpolation.verifySettings()

            if type(self.ripplemode) == bool:
                self.setRippleMode(self.ripplemode)
            elif type(self.ripplemode) != int:
                raise EquationException("{}: Invalid type of ripple mode option: {}".format(self.name, type(self.ripplemode)))
            else:
                opt = [RIPPLE_MODE_NEGLECT, RIPPLE_MODE_BOX, RIPPLE_MODE_GAUSSIAN]
                if self.ripplemode not in opt:
                    raise EquationException("{}: Invalid option for ripple mode: {}.".format(self.name, self.ripplemode))
 
            if type(self.synchrotronmode) == bool:
                self.setSynchrotronMode(self.synchrotronmode)
            elif type(self.synchrotronmode) != int:
                raise EquationException("{}: Invalid type of synchrotron mode option: {}".format(self.name, type(self.synchrotronmode)))
            else:
                opt = [SYNCHROTRON_MODE_NEGLECT, SYNCHROTRON_MODE_INCLUDE]
                if self.synchrotronmode not in opt:
                    raise EquationException("{}: Invalid option for synchrotron mode: {}".format(self.name, self.synchrotronmode))

            if type(self.timevaryingbmode) == bool:
                self.setTimeVaryingBMode(self.timevaryingbmode)
            elif type(self.timevaryingbmode) != int:
                raise EquationException(f"{self.name}: Invalid type of time-varying B mode option: {self.timevaryingbmode}.")
            else:
                opt = [TIME_VARYING_B_MODE_NEGLECT, TIME_VARYING_B_MODE_INCLUDE]
                if self.timevaryingbmode not in opt:
                    raise EquationException(f"{self.name}: Invalid option for time-varying B mode: {self.timevaryingbmode}.")

            self.transport.verifySettings()
        elif self.mode != DISTRIBUTION_MODE_NUMERICAL:
            # if fluid mode and analytical distribution,
            # initial profiles must be provided:
            self.verifyInitialProfiles()


    def verifyInitialDistribution(self):
        """
        Verifies that the initial distribution function has
        been set correctly and consistently.
        """
        if self.init is None:
            raise EquationException("{}: No initial distribution function specified.".format(self.name))

        nr = self.init['r'].size
        p1, p2 = None, None
        p1name, p2name = None, None
        np1, np2 = 0, 0

        if self.init['p'].size > 0 and self.init['xi'].size > 0:
            p1name = 'p'
            p2name = 'xi'
        elif self.init['ppar'].size > 0 and self.init['pperp'].size > 0:
            p1name = 'ppar'
            p2name = 'pperp'
        else:
            raise EquationException("{}: No momentum grid given for initial value.".format(self.name))

        p1 = self.init[p1name]
        p2 = self.init[p2name]

        if len(p1.shape) != 1:
            raise EquationException("{}: Invalid dimensions of momentum grid '{}'. Must be 1D array.".format(self.name, p1name))
        elif len(p2.shape) != 1:
            raise EquationException("{}: Invalid dimensions of momentum grid '{}'. Must be 1D array.".format(self.name, p2name))

        np1 = p1.size
        np2 = p2.size

        if self.init['x'].shape != (nr, np2, np1):
            raise EquationException("{}: Invalid size of initial distribution function: {}. Expected: {}.".format(self.name, self.init['x'].shape, (nr, np2, np1)))


    def verifyInitialProfiles(self):
        """
        Verifies that the initial density and temperature profiles
        are set correctly.
        """
        if (self.n0 is None) or (self.T0 is None):
            raise EquationException("{}: No initial density and/or temperature profiles specified.".format(self.name))
        if (self.rn0 is None) or (self.rT0 is None):
            raise EquationException("{}: No radial grids specified for the density and/or temperature profiles.".format(self.name))

        if (self.n0.ndim != 1) or (self.rn0.ndim != 1) or (self.n0.size != self.rn0.size):
            raise EquationException("{}: Invalid number of elements of density profile: {}. Corresponding radial grid has {} elements."
                .format(self.name, self.n0.size, self.rn0.size))
        if (self.T0.ndim != 1) or (self.rT0.ndim != 1) or (self.T0.size != self.rT0.size):
            raise EquationException("{}: Invalid number of elements of temperature profile: {}. Corresponding radial grid has {} elements."
                .format(self.name, self.T0.size, self.rT0.size))


