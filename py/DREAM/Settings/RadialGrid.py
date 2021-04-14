#
# RadialGrid settings object.
########################################

import numpy as np
import matplotlib.pyplot as plt
import pathlib
import scipy.interpolate
from DREAM.DREAMException import DREAMException
from .Equations.EquationException import EquationException
from .LUKEMagneticField import LUKEMagneticField


TYPE_CYLINDRICAL = 1
TYPE_ANALYTIC_TOROIDAL = 2
TYPE_NUMERICAL = 3

# Numerical magnetic field file formats
FILE_FORMAT_LUKE = 1


class RadialGrid:
    

    def __init__(self, ttype=1):
        """
        Constructor.
        """
        self.setType(ttype=ttype)

        # Cylindrical settings
        self.a  = 0.0
        self.b  = 0.0
        self.B0 = 0.0
        self.nr = int(0)
        self.r0 = 0.0

        # Analytic toroidal settings
        self.R0 = 2.0
        self.ntheta = 20
        self.Delta = None       # Shafranov shift
        self.Delta_r = None
        self.delta = None       # Triangularity
        self.delta_r = None
        self.GOverR0 = None     # R*Bphi/R0
        self.G_r = None
        self.kappa = None       # Elongation
        self.kappa_r = None
        self.psi_p0 = None      # Reference poloidal flux
        self.psi_p0_r = None
        # Ripple parameters
        self.ripple_ncoils = 0
        self.ripple_deltacoils = 0.0
        self.ripple_m = None
        self.ripple_n = None
        self.ripple_dB_B = None
        self.ripple_r = None
        self.ripple_t = None

        # Numerical magnetic field parameters
        self.num_filename = None
        self.num_fileformat = None
        self.num_magneticfield = None   # Magnetic field class parsing data

        # prescribed arbitrary grid
        self.r_f = None 


    #######################
    # SETTERS
    #######################
    def setCustomGridPoints(self, r_f):
        """
        (Cylindrical, Analytic toroidal)
        Set an arbitrary custom grid point distribution
        on the radial flux grid (i.e. the locations of
        the cell edges). This overrides the grid resolution
        'nr', which will be taken as the number of cells
        described by the prescribed grid points.

        :param float r_f: List of radial flux grid points
        """
        if self.nr != 0 or self.a != 0 or self.r0 != 0:
            #raise EquationException("RadialGrid: Cannot assign custom grid points while prescribing 'nr', 'a' or 'r0'.")         
            print("*WARNING* RadialGrid: Prescibing custom radial grid overrides 'nr', 'a' and 'r0'.")
            self.nr = int(0)
            self.a  = 0.0
            self.r0 = 0.0

        if type(r_f)==list:
            r_f = np.array(r_f)
        if np.size(r_f)<2:
            raise EquationException("RadialGrid: Custom grid point vector 'r_f' must have size 2 or greater.")
        for i in range(np.size(r_f)-1):
            if r_f[i+1]<r_f[i]:
                raise EquationException("RadialGrid: Custom grid points 'r_f' must be an array of increasing numbers.")
        if np.min(r_f)<0:
            raise EquationException("RadialGrid: Custom grid points must be non-negative.")
        self.r_f = r_f

    def setB0(self, B0):
        """
        (Cylindrical)
        Set the on-axis magnetic field strength.
        """
        if B0 <= 0:
            raise EquationException("RadialGrid: Invalid value assigned to 'B0': {}. Must be >0.".format(B0))
        
        self.B0 = float(B0)


    def setInnerRadius(self, r0):
        """
        (Cylindrical, Analytic toroidal)
        Set the innermost radial point to simulate.
        """
        if r0 < 0:
            raise EquationException("RadialGrid: Invalid value assigned to innermost radius 'r0': {}".format(r0))
        if self.r_f is not None:
            print("*WARNING* RadialGrid: Prescibing 'Inner radius' r0 overrides the custom radial grid 'r_f'.")
            self.r_f = None

        self.r0 = r0


    def setMinorRadius(self, a):
        """
        (Cylindrical, Analytic toroidal)
        Set the plasma minor radius.
        """
        if a <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(a))
        if self.r_f is not None:
            print("*WARNING* RadialGrid: Prescibing 'Minor radius' a overrides the custom radial grid 'r_f'.")
            self.r_f = None

        self.a = float(a)


    def setMajorRadius(self, R0):
        """
        (Analytic toroidal)
        Set the tokamak major radius.
        """
        if R0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to major radius 'R0': {}".format(R0))

        self.R0 = float(R0)

    def setWallRadius(self, wall_radius):
        """
        (Cylindrical, Analytic toroidal)
        Set the minor radius of the wall
        """
        self.b = float(wall_radius)

    def setNr(self, nr):
        """
        (Cylindrical, Analytic toroidal)
        Set the number of radial grid points to use.
        """
        if nr <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to 'nr': {}".format(nr))
        if self.r_f is not None:
            print("*WARNING* RadialGrid: Prescibing 'Nr' overrides the custom radial grid 'r_f'.")
            self.r_f = None
            
        self.nr = int(nr)


    def setNtheta(self, ntheta):
        """
        (Analytic toroidal and numerical)
        Set the number of grid points to use for the poloidal grid on which bounce
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
        if name not in ['Delta', 'delta', 'GOverR0', 'kappa', 'psi_p0']:
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


    def setShaping(self, psi, GOverR0, rpsi=0.0, rG=0.0,
        Delta=0.0, rDelta=0.0, delta=0.0, rdelta=0.0,
        kappa=1.0, rkappa=0.0):
        """
        (Analytic toroidal)
        Set the plasma shape parameters to use with the magnetic field.

        :param GOverR0:      Toroidal magnetic field component, ``R*Bphi``, normalized by ``R0``.
        :param rG:     Radial grid for ``GOverR0``.
        :param psi:    Reference poloidal flux, normalized by ``R0``.
        :param rpsi:   Radial grid for ``psi``.
        :param Delta:  Shafranov shift.
        :param rDelta: Radial grid for Shafranov shift.
        :param delta:  Triangularity.
        :param rdelta: Radial grid for triangularity.
        :param kappa:  Elongation.
        :param rkappa: Radial grid for elongation.
        """
        self.setShapeParameter('Delta',   r=rDelta, data=Delta)
        self.setShapeParameter('delta',   r=rdelta, data=delta)
        self.setShapeParameter('GOverR0', r=rG,     data=GOverR0)
        self.setShapeParameter('kappa',   r=rkappa, data=kappa)
        self.setShapeParameter('psi_p0',  r=rpsi,   data=psi)


    def setRipple(self, m, n, dB_B, ncoils=0, deltacoils=0, r=[0], t=[0]):
        """
        Enable the ripple pitch scattering term.

        :param list m:           Poloidal mode numbers of magnetic perturbation(s).
        :param list n:           Toroidal mode numbers of magnetic perturbation(s).
        :param dB_B:             Magnetic perturbations (shape: nModes x nt x nr).
        :param int ncoils:       Number of toroidal field coils.
        :param float deltacoils: Distance between toroidal field coils.
        :param r:                Radial grid on which the magnetic perturbations are given.
        :param t:                Time grid on which the magnetic perturbations are given.
        """
        if type(m) == list: m = np.array(m)
        elif np.isscalar(m): m = np.array([float(m)])

        if type(n) == list: n = np.array(n)
        elif np.isscalar(n): n = np.array([float(n)])

        if type(r) == list: r = np.array(r)
        elif np.isscalar(r): r = np.array([float(r)])

        if type(t) == list: t = np.array(t)
        elif np.isscalar(t): t = np.array([float(t)])

        if type(dB_B) == list:
            dB_B = np.array(dB_B)
        if type(dB_B) == float or dB_B.ndim == 1:
            dB_B = np.ones((m.size, t.size, r.size)) * dB_B

        if m.size != n.size:
            raise EquationException("RadialGrid: m and n must have the same number of elements.")
        elif dB_B.ndim == 1 and dB_B.size == m.size:
            dB_B = dB_B*np.ones((m.size, t.size, r.size))
        elif dB_B.ndim != 3 or dB_B.shape != (m.size, t.size, r.size):
            raise EquationException("RadialGrid: Invalid dimensions of parameter 'dB_B'. Expected {}, but array has {}.".format((m.size, t.size, r.size), dB_B.shape))

        self.ripple_ncoils = int(ncoils)
        self.ripple_deltacoils = float(deltacoils)
        self.ripple_m = m
        self.ripple_n = n
        self.ripple_dB_B = dB_B
        self.ripple_r = r
        self.ripple_t = t


    def setNumerical(self, filename, format=FILE_FORMAT_LUKE):
        """
        Sets the numerical magnetic field to use for the simulation.

        :param str filename: Name of file containing magnetic field data.
        :param int format:   Format of the magnetic field data in the given file.
        """
        self.type = TYPE_NUMERICAL
        self.num_filename = filename
        
        if format is not None:
            self.num_fileformat = format

        if format == FILE_FORMAT_LUKE:
            self.num_magneticfield = LUKEMagneticField(filename)

        self.a = self.num_magneticfield.a


    def setType(self, ttype):
        """
        Set the type of radial grid to use.
        """
        types = [TYPE_CYLINDRICAL, TYPE_ANALYTIC_TOROIDAL, TYPE_NUMERICAL]
        if ttype in types:
            self.type = ttype
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(ttype))


    def visualize(self, *args, ax=None, show=None, **kwargs):
        """
        Visualize the current magnetic field.

        :param int nr:     Number of flux surfaces to show.
        :param int ntheta: Number of poloidal angles per flux surface.
        """
        # Ensure that settings are valid...
        self.verifySettings()

        if self.type == TYPE_ANALYTIC_TOROIDAL:
            self.visualize_analytic(*args, ax=ax, show=show, **kwargs)
        elif self.type == TYPE_NUMERICAL:
            self.num_magneticfield.visualize(*args, ax=ax, show=show, **kwargs)
        else:
            raise DREAMException("RadialGrid: Can only visualize the analytic toroidal magnetic field.")

        
    def visualize_analytic(self, nr=10, ntheta=40, ax=None, show=None):
        """
        Visualize an analytic toroidal magnetic field.
        """
        red   = (249/255, 65/255, 68/255)
        black = (87/255, 117/255, 144/255)
        gray  = (190/255, 190/255, 190/255)

        # Set up axes (if not already done)
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        # Generate r/theta grid
        if self.a == 0:
            # Custom radial grid (i.e. not generated by specifying r0 and a)
            r_f = self.r_f
            # Resample
            r_f = np.interp(np.linspace(0, 1, nr)*r_f.size, range(r_f.size), r_f)
        else:
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
                f = scipy.interpolate.interp1d(rParam, param, kind='cubic', bounds_error=False, fill_value='extrapolate')

            return f(r)

        # Interpolate shaping parameters
        Delta = lambda r : interppar(r, self.Delta_r, self.Delta)
        delta = lambda r : interppar(r, self.delta_r, self.delta)
        kappa = lambda r : interppar(r, self.kappa_r, self.kappa)

        # Construct flux surfaces
        R = lambda r, t : self.R0 + Delta(r) + r*np.cos(t + delta(r)*np.sin(t))
        Z = lambda r, t : r*kappa(r)*np.sin(t)

        # Flux surfaces
        ax.plot(R(rr, tt), Z(rr, tt), color=gray, linewidth=1)
        # Limiter
        ax.plot(R(r_f[-1], tt), Z(r_f[-1], tt), color=black, linewidth=2)
        # Wall
        ax.plot(R(np.array([self.b]), tt), Z(np.array([self.b]), tt), color=red, linewidth=2)
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
        def scal(v):
            if type(v) == np.ndarray: return v[0]
            else: return v

        self.type = data['type']

        if 'wall_radius' in data:
            self.b = data['wall_radius']
            if type(self.b) == np.ndarray:
                self.b = float(self.b[0])
            else:
                self.b = float(self.b)

        if self.type == TYPE_CYLINDRICAL or self.type == TYPE_ANALYTIC_TOROIDAL or self.type == TYPE_NUMERICAL:
            self.a = data['a']
            self.nr = data['nr']
            self.r0 = data['r0']
            if 'r_f' in data:
                self.r_f = data['r_f']

        if self.type == TYPE_CYLINDRICAL:
            self.B0 = data['B0']
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            self.R0 = data['R0']
            self.ntheta = data['ntheta']

            self.Delta = data['Delta']['x']
            self.Delta_r = data['Delta']['r']
            self.delta = data['delta']['x']
            self.delta_r = data['delta']['r']
            self.GOverR0 = data['GOverR0']['x']
            self.G_r = data['GOverR0']['r']
            self.kappa = data['kappa']['x']
            self.kappa_r = data['kappa']['r']
            self.psi_p0 = data['psi_p0']['x']
            self.psi_p0_r = data['psi_p0']['r']
        elif self.type == TYPE_NUMERICAL:
            self.num_filename = data['filename']
            self.ntheta = data['ntheta']

            if 'fileformat' in data:
                self.num_fileformat = data['fileformat']
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        if 'ripple' in data:
            self.ripple_ncoils = int(scal(data['ripple']['ncoils']))
            self.ripple_deltacoils = float(scal(data['ripple']['deltacoils']))
            self.ripple_m = data['ripple']['m']
            self.ripple_n = data['ripple']['n']
            self.ripple_dB_B = data['ripple']['x']
            self.ripple_r = data['ripple']['r']
            self.ripple_t = data['ripple']['t']


    def todict(self, verify=True):
        """
        Returns the settings in this object as a Python dictionary.
        """
        if verify:
            self.verifySettings()

        data = {
            'type': self.type
        }

        if self.type == TYPE_CYLINDRICAL or self.type == TYPE_ANALYTIC_TOROIDAL or self.type == TYPE_NUMERICAL:
            data['a'] = self.a
            data['nr'] = self.nr
            data['r0'] = self.r0
            data['wall_radius'] = self.b
            if self.r_f is not None:
                data['r_f'] = self.r_f

        if self.type == TYPE_CYLINDRICAL:
            data['B0'] = self.B0
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            data['R0'] = self.R0
            data['ntheta'] = self.ntheta

            data['Delta']   = {'x': self.Delta, 'r': self.Delta_r}
            data['delta']   = {'x': self.delta, 'r': self.delta_r}
            data['GOverR0'] = {'x': self.GOverR0, 'r': self.G_r}
            data['kappa']   = {'x': self.kappa, 'r': self.kappa_r}
            data['psi_p0']  = {'x': self.psi_p0, 'r': self.psi_p0_r}
        elif self.type == TYPE_NUMERICAL:
            data['filename'] = self.num_filename
            data['ntheta'] = self.ntheta

            if self.num_fileformat is not None:
                data['fileformat'] = self.num_fileformat
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        if self.ripple_ncoils > 0 or self.ripple_deltacoils > 0:
            data['ripple'] = {
                'ncoils': self.ripple_ncoils,
                'deltacoils': self.ripple_deltacoils,
                'm': self.ripple_m,
                'n': self.ripple_n,
                'x': self.ripple_dB_B,
                'r': self.ripple_r,
                't': self.ripple_t
            }

        return data
        
            
    def verifySettings(self):
        """
        Verfiy that the RadialGrid settings are consistent.
        """
        types = [TYPE_CYLINDRICAL, TYPE_ANALYTIC_TOROIDAL, TYPE_NUMERICAL]
        if self.type in types:
            if (self.a is None or self.a <= 0) and self.r_f is None:
                raise DREAMException("RadialGrid: Invalid value assigned to minor radius 'a': {}".format(self.a))
            elif (self.r0 is None or self.r0 < 0) and self.r_f is None:
                raise DREAMException("RadialGrid: Invalid value assigned to innermost simulated radius 'r0': {}".format(self.r0))
            elif self.b is None or self.b<self.a:
                raise DREAMException("RadialGrid: Invalid value assigned to wall radius 'b' (must be explicitly set to >= 'a' using 'setWallRadius'): ".format(self.b))
            if self.r0 >= self.a and self.r_f is None:
                raise DREAMException("RadialGrid: 'r0' must be strictly less than 'a'.")
            if self.nr <= 0 and self.r_f is None:
                raise DREAMException("RadialGrid: Invalid value assigned 'nr': {}. Must be > 0.".format(self.nr))
            if not np.isscalar(self.b):
                raise DREAMException("RadialGrid: The specified wall radius is not a scalar: {}.".format(self.b))

        if self.type == TYPE_CYLINDRICAL:
            if self.B0 is None or self.B0 <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'B0': {}".format(self.B0))
        elif self.type == TYPE_ANALYTIC_TOROIDAL:
            if self.R0 is None or self.R0 <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to tokamak major radius 'R0': {}".format(self.R0))
            elif self.ntheta <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'ntheta': {}. Must be > 0.".format(self.ntheta))

            self.verifySettingsShapeParameter('Delta')
            self.verifySettingsShapeParameter('delta')
            self.verifySettingsShapeParameter('GOverR0')
            self.verifySettingsShapeParameter('kappa')
            self.verifySettingsShapeParameter('psi_p0')

            if np.size(self.Delta_r)>1:
                if self.Delta_r[0]==0 and self.Delta[0]!=0:
                    print("*WARNING* RadialGrid: Shape parameter 'Delta' (Shafranov shift) is non-zero at r=0, which is inconsistent (add Delta(0) to the major radius R0 instead)")
            elif self.Delta!=0:
                print("*WARNING* RadialGrid: Shape parameter 'Delta' (Shafranov shift) is assigned a constant non-zero value. It is recommended to add its value to the major radius R0 instead")
            if np.size(self.delta_r)>1:
                if self.delta_r[0]==0 and self.delta[0]!=0:
                    print("*WARNING* RadialGrid: Shape parameter 'delta' (triangularity) is non-zero at r=0, which is inconsistent with Grad-Shafranov")
        elif self.type == TYPE_NUMERICAL:
            if type(self.num_filename) != str:
                raise DREAMException("RadialGrid: No numerical magnetic field file specified.")
            elif not pathlib.Path(self.num_filename).is_file():
                raise DREAMException("RadialGrid: The specified numerical magnetic field file does not exist.")
            elif self.ntheta <= 0:
                raise DREAMException("RadialGrid: Invalid value assigned to 'ntheta': {}. Must be > 0.".format(self.ntheta))

            formats = [FILE_FORMAT_LUKE]
            if (self.num_fileformat is not None) and (self.num_fileformat not in formats):
                raise DREAMException("RadialGrid: Unrecognized file format specified for numerical magnetic field: {}.".format(self.num_fileformat))
        else:
            raise DREAMException("RadialGrid: Unrecognized grid type specified: {}.".format(self.type))

        # Ripple settings
        if self.ripple_ncoils > 0 or self.ripple_deltacoils > 0:
            if type(self.ripple_m) != np.ndarray or self.ripple_m.ndim != 1:
                raise EquationException("RadialGrid: Invalid type or shape of 'ripple_m'.")
            elif type(self.ripple_n) != np.ndarray or self.ripple_n.ndim != 1:
                raise EquationException("RadialGrid: Invalid type or shape of 'ripple_n'.")
            elif self.ripple_m.size != self.ripple_n.size:
                raise EquationException("RadialGrid: 'ripple_m' and 'ripple_n' must have the same number of elements.")
            elif type(self.ripple_r) != np.ndarray or self.ripple_r.ndim != 1:
                raise EquationException("RadialGrid: Invalid type or shape of 'ripple_r'.")
            elif type(self.ripple_t) != np.ndarray or self.ripple_t.ndim != 1:
                raise EquationException("RadialGrid: Invalid type or shape of 'ripple_t'.")
            elif type(self.ripple_dB_B) != np.ndarray or self.ripple_dB_B.shape != (self.ripple_m.size, self.ripple_r.size, self.ripple_t.size):
                raise EquationException("RadialGrid: Invalid type or shape of 'ripple_dB_B'.".format(self.ripple_dB_B))
        
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


