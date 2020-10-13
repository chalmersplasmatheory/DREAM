#
# RadialGrid settings object.
########################################

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from DREAM.DREAMException import DREAMException


TYPE_CYLINDRICAL = 1
TYPE_ANALYTIC_TOROIDAL = 2


class RadialGrid:
    

    def __init__(self, ttype=1):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Cylindrical settings
        self.a  = 0.0
        self.B0 = 0.0
        self.nr = int(0)
        self.r0 = 0.0

        # Analytic toroidal settings
        self.R0 = 2.0
        self.ntheta = 30
        self.Delta = None       # Shafranov shift
        self.Delta_r = None
        self.delta = None       # Triangularity
        self.delta_r = None
        self.G = None           # R*Bphi
        self.G_r = None
        self.kappa = None       # Elongation
        self.kappa_r = None
        self.psi_p0 = None      # Reference poloidal flux
        self.psi_p0_r = None


    #######################
    # SETTERS
    #######################
    def setB0(self, B0):
        """
        (Cylindrical)
        Set the on-axis magnetic field strength.
        """
        if B0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'B0'.")
        
        self.B0 = float(B0)


    def setInnerRadius(self, r0):
        """
        (Cylindrical, Analytic toroidal)
        Set the innermost radial point to simulate.
        """
        if r0 < 0:
            raise DREAMException("RadialGrid: Invalid value assigned to innermost radius 'r0': {}".format(r0))

        self.r0 = r0


    def setMinorRadius(self, a):
        """
        (Cylindrical, Analytic toroidal)
        Set the plasma minor radius.
        """
        if a <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(a))

        self.a = float(a)


    def setMajorRadius(self, R0):
        """
        (Analytic toroidal)
        Set the tokamak major radius.
        """
        if R0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to major radius 'R0': {}".format(R0))

        self.R0 = float(R0)


    def setNr(self, nr):
        """
        (Cylindrical, Analytic toroidal)
        Set the number of radial grid points to use.
        """
        if nr <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'nr': {}".format(nr))

        self.nr = int(nr)


    def setNtheta(self, ntheta):
        """
        (Analytic toroidal)
        Set the number of grid points to use on the poloidal on which bounce
        averages are calculated.
        """
        if ntheta <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'ntheta': {}".format(ntheta))

        self.ntheta = ntheta


    def setShapeParameter(self, name, data, r=0.0):
        """
        (Analytic toroidal)
        Set a specific magnetic field shape parameter.
        """
        if name not in ['Delta', 'delta', 'G', 'kappa', 'psi_p0']:
            raise DREAMException("RadialGrid: Invalid name of shape parameter specified: '{}'.".format(name))

        if type(data) in [float, int]:
            data = np.array([float(data)])
        elif type(data) == list:
            data = np.array(data)

        if type(r) in [float, int]:
            r = np.array([float(r)])
        elif type(r) == list:
            r = np.array(r)

        setattr(self, name, data)
        setattr(self, name+'_r', r)


    def setShaping(self, psi, G, rpsi=0.0, rG=0.0,
        Delta=0.0, rDelta=0.0, delta=0.0, rdelta=0.0,
        kappa=1.0, rkappa=0.0):
        """
        (Analytic toroidal)
        Set the plasma shape parameters to use with the magnetic field.

        :param G:      Toroidal magnetic field component, R*Bphi.
        :param rG:     Radial grid for ``G``.
        :param psi:    Reference poloidal flux.
        :param rpsi:   Radial grid for ``psi``.
        :param Delta:  Shafranov shift.
        :param rDelta: Radial grid for Shafranov shift.
        :param delta:  Triangularity.
        :param rdelta: Radial grid for triangularity.
        :param kappa:  Elongation.
        :param rkappa: Radial grid for elongation.
        """
        self.setShapeParameter('Delta',  r=rDelta, data=Delta)
        self.setShapeParameter('delta',  r=rdelta, data=delta)
        self.setShapeParameter('G',      r=rG,     data=G)
        self.setShapeParameter('kappa',  r=rkappa, data=kappa)
        self.setShapeParameter('psi_p0', r=rpsi,   data=psi)


    def setType(self, ttype):
        """
        Set the type of radial grid to use.
        """
        if ttype == TYPE_CYLINDRICAL:
            self.type = ttype
        elif ttype == TYPE_ANALYTIC_TOROIDAL:
            self.type = ttype
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(ttype))


    def visualize(self, nr=10, ntheta=30, ax=None, show=None):
        """
        Visualize the current magnetic field.

        :param int nr:     Number of flux surfaces to show.
        :param int ntheta: Number of poloidal angles per flux surface.
        """
        if self.type != TYPE_ANALYTIC_TOROIDAL:
            raise DREAMException("RadialGrid: Can only visualize the analytic toroidal magnetic field.")

        # Ensure that settings are valid...
        self.verifySettings()

        # Set up axes (if not already done)
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        # Generate r/theta grid
        r_f = np.linspace(self.r0, self.a, nr+1)
        r   = (r_f[:-1] + r_f[1:]) / 2.0
        t   = np.linspace(0, 2*np.pi, ntheta)

        rr, tt = np.meshgrid(r, t)

        # Helper interpolation function
        def interppar(r, rParam, param):
            f = None
            if param.size == 1:
                return param[0] * np.ones(r.shape)
            elif param.size == 2:
                f = scipy.interpolate.interp1d(rParam, param)
            else:
                # Not exactly what is done in the kernel (which uses Steffen
                # interpolation, which is cubic and guarantees positivity)
                f = scipy.interpolate.interp1d(rParam, param, kind='cubic')

            return f(r)

        # Interpolate shaping parameters
        Delta = lambda r : interppar(r, self.Delta_r, self.Delta)
        delta = lambda r : interppar(r, self.delta_r, self.delta)
        kappa = lambda r : interppar(r, self.kappa_r, self.kappa)

        # Construct flux surfaces
        R = lambda r, t : self.R0 + Delta(r) + r*np.cos(t + delta(r)*np.sin(t))
        Z = lambda r, t : r*kappa(r)*np.sin(t)

        # Flux surfaces
        ax.plot(R(rr, tt), Z(rr, tt), 'k', linewidth=2)
        # Limiter
        ax.plot(R(r_f[-1], tt), Z(r_f[-1], tt), 'r', linewidth=3)
        ax.axis('equal')

        ax.set_xlabel('Major radius $R$ (m)')
        ax.set_ylabel('Height $Z$ (m)')

        if show:
            plt.show()

        return ax

    
    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.type = data['type']

        if self.type == TYPE_CYLINDRICAL:
            self.a = data['a']
            self.B0 = data['B0']
            self.nr = data['nr']
            self.r0 = data['r0']
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            self.a  = data['a']
            self.nr = data['nr']
            self.r0 = data['r0']
            self.R0 = data['R0']
            self.ntheta = data['ntheta']

            self.Delta = data['Delta']['x']
            self.Delta_r = data['Delta']['r']
            self.delta = data['delta']['x']
            self.delta_r = data['delta']['r']
            self.G = data['G']['x']
            self.G_r = data['G']['r']
            self.kappa = data['kappa']['x']
            self.kappa_r = data['kappa']['r']
            self.psi_p0 = data['psi_p0']['x']
            self.psi_p0_r = data['psi_p0']['r']
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))


    def todict(self, verify=True):
        """
        Returns the settings in this object as a Python dictionary.
        """
        if verify:
            self.verifySettings()

        data = {
            'type': self.type
        }

        if self.type == TYPE_CYLINDRICAL:
            data['a'] = self.a
            data['B0'] = self.B0
            data['nr'] = self.nr
            data['r0'] = self.r0
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            data['a'] = self.a
            data['nr'] = self.nr
            data['r0'] = self.r0
            data['R0'] = self.R0
            data['ntheta'] = self.ntheta

            data['Delta']  = {'x': self.Delta, 'r': self.Delta_r}
            data['delta']  = {'x': self.delta, 'r': self.delta_r}
            data['G']      = {'x': self.G, 'r': self.G_r}
            data['kappa']  = {'x': self.kappa, 'r': self.kappa_r}
            data['psi_p0'] = {'x': self.psi_p0, 'r': self.psi_p0_r}
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        return data
        
            
    def verifySettings(self):
        """
        Verfiy that the RadialGrid settings are consistent.
        """
        if self.type == TYPE_CYLINDRICAL:
            if self.a is None or self.a <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(self.a))
            elif self.B0 is None or self.B0 <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'B0': {}".format(self.B0))
            elif self.r0 is None or self.r0 < 0:
                raise DREAMException("RadialGrid: Invalid value assigned to innermost simulated radius 'r0': {}".format(self.r0))

            if self.nr <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned 'nr': {}. Must be > 0.".format(self.nr))

            if self.r0 >= self.a:
                raise DREAMException("RadialGrid: 'r0' must be strictly less than 'a'.")
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            if self.a is None or self.a <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(self.a))
            elif self.r0 is None or self.r0 < 0:
                raise DREAMException("RadialGrid: Invalid value assigned to innermost simulated radius 'r0': {}".format(self.r0))
            elif self.R0 is None or self.R0 <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to tokamak major radius 'R0': {}".format(self.R0))

            if self.nr <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned 'nr': {}. Must be > 0.".format(self.nr))
            elif self.ntheta <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'ntheta': {}. Must be > 0.".format(self.ntheta))

            if self.r0 >= self.a:
                raise DREAMException("RadialGrid: 'r0' must be strictly less than 'a'.")

            self.verifySettingsShapeParameter('Delta')
            self.verifySettingsShapeParameter('delta')
            self.verifySettingsShapeParameter('G')
            self.verifySettingsShapeParameter('kappa')
            self.verifySettingsShapeParameter('psi_p0')
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))


    def verifySettingsShapeParameter(self, shapeparam):
        """
        Verify the settings of the named shape parameter.
        
        :param str shapeparam: Name of shape parameter to verify settings for.
        """
        v = getattr(self, shapeparam)
        r = getattr(self, shapeparam+'_r')

        if v is None or type(v) != np.ndarray:
            raise DREAMException("RadialGrid: Invalid type of shape parameter '{}': {}.".format(shapeparam, type(v)))
        elif r is None or type(r) != np.ndarray:
            raise DREAMException("RadialGrid: Invalid type of radial grid for shape parameter '{}': {}.".format(shapeparam, type(r)))

        if v.shape != r.shape:
            raise DREAMException("RadialGrid: Dimensions mismatch between shape parameter '{}' {} and its radial grid {}.".format(shapeparam, v.shape, r.shape))
        

