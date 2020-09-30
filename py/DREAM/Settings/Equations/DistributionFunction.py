
import numpy as np
from DREAM.Settings.Equations.EquationException import EquationException
from . UnknownQuantity import UnknownQuantity
from .. TransportSettings import TransportSettings


# BOUNDARY CONDITIONS (WHEN f_re IS DISABLED)
BC_F_0        = 1
BC_PHI_CONST  = 2
BC_DPHI_CONST = 3

# Interpolation methods for advection term in kinetic equation
AD_INTERP_CENTRED  = 1
AD_INTERP_UPWIND   = 2
AD_INTERP_UPWIND_2ND_ORDER = 3
AD_INTERP_DOWNWIND = 4
AD_INTERP_QUICK    = 5
AD_INTERP_SMART    = 6
AD_INTERP_MUSCL    = 7
AD_INTERP_OSPRE    = 8
AD_INTERP_TCDF     = 9


class DistributionFunction(UnknownQuantity):
    

    def __init__(self, settings, name, grid,
        f=None, initr=None, initp=None, initxi=None,
        initppar=None, initpperp=None,
        rn0=None, n0=None, rT0=None, T0=None, bc=BC_PHI_CONST,
        ad_int_r=AD_INTERP_CENTRED, ad_int_p1=AD_INTERP_CENTRED,
        ad_int_p2=AD_INTERP_CENTRED, fluxlimiterdamping=1.0):
        """
        Constructor.
        """
        super().__init__(settings=settings)

        self.name = name
        self.grid = grid

        self.boundarycondition = bc
        
        self.transport = TransportSettings(kinetic=True)

        self.adv_interp_r  = ad_int_r 
        self.adv_interp_p1 = ad_int_p1
        self.adv_interp_p2 = ad_int_p2 
        self.fluxlimiterdamping = fluxlimiterdamping

        self.n0  = rn0
        self.rn0 = n0

        self.T0  = rT0
        self.rT0 = T0

        self.init = None

        if f is not None:
            self.setInitialValue(f, r=initr, p=initp, xi=initxi, ppar=initppar, pperp=initpperp)
        elif n0 is not None:
            self.setInitialProfiles(rn0=rn0, n0=n0, rT0=rT0, T0=T0)


    def setBoundaryCondition(self, bc):
        """
        Sets the boundary condition at p=pmax to use when 'f_re' is disabled.
        When 'f_re' is enabled, only one boundary condition can be used, in
        which case this flag is ignored.
        """
        self.boundarycondition = bc


    def setAdvectionInterpolationMethod(self,ad_int=None, ad_int_r=AD_INTERP_CENTRED,
        ad_int_p1=AD_INTERP_CENTRED,ad_int_p2=AD_INTERP_CENTRED,fluxlimiterdamping=1.0):
        """
        Sets the interpolation method that is used
        in the advection terms of the kinetic equation.
        To set all three components, provide ad_int. Otherwise
        the three components can use separate interpolation methods.
        """
        self.fluxlimiterdamping = fluxlimiterdamping
        if ad_int is not None:
            self.adv_interp_r  = ad_int
            self.adv_interp_p1 = ad_int
            self.adv_interp_p2 = ad_int
        else:
            self.adv_interp_r  = ad_int_r
            self.adv_interp_p1 = ad_int_p1
            self.adv_interp_p2 = ad_int_p2


    def setInitialProfiles(self, n0, T0, rn0=None, rT0=None):
        """
        Sets the initial density and temperature profiles of the
        electron population.

        rn0: Radial grid on which the density is given.
        n0:  Electron density profile.
        rT0: Radial grid on which the temperature is given.
        T0:  Electron temperature profile.
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
        Set the initial value of this electron distribution function.

        f:     Array representing the distribution function value on the grid
               (must have size (nr, nxi, np) or (nr, npperp, nppar))
        r:     Radial grid on which the initial distribution is given.
          MOMENTUM GRID SETTINGS:
        p:     Momentum grid.
        xi:    Pitch grid.
          OR
        ppar:  Parallel momentum grid.
        pperp: Perpendicular momentum grid.
        """
        self.init = {}

        if p is not None and xi is not None:
            self.init['x'] = np.asarray(f)
            self.init['r'] = np.asarray(r)
            self.init['p'] = np.asarray(p)
            self.init['xi'] = np.asarray(xi)
            self.init['ppar'] = np.array([])
            self.init['pperp'] = np.array([])
        elif ppar is not None and pperp is not None:
            self.init['x'] = np.asarray(f)
            self.init['r'] = np.asarray(r)
            self.init['ppar'] = np.asarray(ppar)
            self.init['pperp'] = np.asarray(pperp)
            self.init['p'] = np.array([])
            self.init['xi'] = np.array([])
        else:
            raise EquationException("{}: No momentum grid given for initial value.".format(self.name))

        # Reset initial profiles (if any)
        self.rn0 = self.rT0 = None
        self.n0 = self.T0 = None

        self.verifyInitialDistribution()


    def fromdict(self, data):
        """
        Load data for this object from the given dictionary.
        """
        if 'boundarycondition' in data:
            self.boundarycondition = data['boundarycondition']
        if 'adv_interp' in data:
            self.adv_interp_r = data['adv_interp']['r']
            self.adv_interp_p1 = data['adv_interp']['p1']
            self.adv_interp_p2 = data['adv_interp']['p2']
            self.fluxlimiterdamping = data['adv_interp']['fluxlimiterdamping']
        if 'init' in data:
            self.init = data['init']
        elif ('n0' in data) and ('T0' in data):
            self.rn0 = data['n0']['r']
            self.n0  = data['n0']['x']
            self.rT0 = data['T0']['r']
            self.T0  = data['T0']['x']
        else:
            raise EquationException("{}: Unrecognized specification of initial distribution function.".format(self.name))

        if 'transport' in data:
            self.transport.fromdict(data['transport'])
            
        self.verifySettings()


    def todict(self):
        """
        Returns a Python dictionary containing all settings of
        this DistributionFunction object.
        """
        data = {}
        if self.grid.enabled:
            data = {'boundarycondition': self.boundarycondition}
            data['adv_interp'] = {}
            data['adv_interp']['r']  = self.adv_interp_r
            data['adv_interp']['p1'] = self.adv_interp_p1
            data['adv_interp']['p2'] = self.adv_interp_p2
            data['adv_interp']['fluxlimiterdamping'] = self.fluxlimiterdamping
            if self.init is not None:
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
            
            data['transport'] = self.transport.todict()

        return data


    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.grid.enabled:
            bc = self.boundarycondition
            if (bc != BC_F_0) and (bc != BC_PHI_CONST) and (bc != BC_DPHI_CONST):
                raise EquationException("{}: Invalid external boundary condition set: {}.".format(self.name, bc))
            ad_int_r = self.adv_interp_r
            if (ad_int_r != AD_INTERP_CENTRED) and (ad_int_r != AD_INTERP_DOWNWIND) and (ad_int_r != AD_INTERP_UPWIND) and (ad_int_r != AD_INTERP_UPWIND_2ND_ORDER) and (ad_int_r != AD_INTERP_QUICK) and (ad_int_r != AD_INTERP_SMART) and (ad_int_r != AD_INTERP_MUSCL)  and (ad_int_r != AD_INTERP_OSPRE) and (ad_int_r != AD_INTERP_TCDF): 
                raise EquationException("{}: Invalid radial interpolation coefficient set: {}.".format(self.name, ad_int_r))
            ad_int_p1 = self.adv_interp_p1
            if (ad_int_p1 != AD_INTERP_CENTRED) and (ad_int_p1 != AD_INTERP_DOWNWIND) and (ad_int_p1 != AD_INTERP_UPWIND) and (ad_int_p1 != AD_INTERP_UPWIND_2ND_ORDER) and (ad_int_p1 != AD_INTERP_QUICK) and (ad_int_p1 != AD_INTERP_SMART) and (ad_int_p1 != AD_INTERP_MUSCL)  and (ad_int_p1 != AD_INTERP_OSPRE) and (ad_int_p1 != AD_INTERP_TCDF):
                raise EquationException("{}: Invalid p1 interpolation coefficient set: {}.".format(self.name, ad_int_p1))
            ad_int_p2 = self.adv_interp_p2
            if (ad_int_p2 != AD_INTERP_CENTRED) and (ad_int_p2 != AD_INTERP_DOWNWIND) and (ad_int_p2 != AD_INTERP_UPWIND) and (ad_int_p2 != AD_INTERP_UPWIND_2ND_ORDER) and (ad_int_p2 != AD_INTERP_QUICK) and (ad_int_p2 != AD_INTERP_SMART) and (ad_int_p2 != AD_INTERP_MUSCL) and (ad_int_p2 != AD_INTERP_OSPRE) and (ad_int_p2 != AD_INTERP_TCDF):
                raise EquationException("{}: Invalid p2 interpolation coefficient set: {}.".format(self.name, ad_int_p2))
            if (self.fluxlimiterdamping<0.0) or (self.fluxlimiterdamping>1.0):
                raise EquationException("{}: Invalid flux limiter damping coefficient: {}. Choose between 0 and 1.".format(self.name, self.fluxlimiterdamping))
            if self.init is not None:
                self.verifyInitialDistribution()
            elif (self.n0 is not None) or (self.T0 is not None):
                self.verifyInitialProfiles()
            else:
                raise EquationException("{}: Invalid/no initial condition set for the distribution function.".format(self.name))

            self.transport.verifySettings()


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

