# 
# Class for loading and working with a magnetic equilibrium stored in a
# GEQDSK file.
#
# Written by: Mathias Hoppe, 2022
#

import h5py
from matplotlib._contour import QuadContourGenerator
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline, RectBivariateSpline
import scipy.optimize


class EqBase:
    

    def __init__(self, filename, override_psilim=False):
        """
        Constructor.
        """
        self.load(filename, override_psilim=override_psilim)


    def _next_value(self, fh):
        """
        Load the next value from the text stream 'fh'.
        """
        pattern = re.compile(r"[ +\-]?\d+(?:\.\d+(?:[Ee][\+\-]\d\d)?)?")

        for line in fh:
            matches = pattern.findall(line)
            for m in matches:
                if "." in m:
                    yield float(m)
                else:
                    yield int(m)


    def cocos(self, nbr=None):
        """
        Returns COCOS factors given a COCOS number.

        Based on "cocos.m" in CHEASEGui, part of the CRPP Toolbox.
        """
        cocos = {}

        if nbr is None:
            nbr = self.cocos_number

        cocos['exp_Bp'] = 1 if nbr >= 11 else 0
        
        if nbr in [1, 11]:
            cocos['sigma_Bp'] = +1
            cocos['sigma_RphiZ'] = +1
            cocos['sigma_rhothetaphi'] = +1
            cocos['sigma_q_pos'] = +1
            cocos['sign_pprime_pos'] = -1
        elif nbr in [2, 12]:
            cocos['sigma_Bp'] = +1
            cocos['sigma_RphiZ'] = -1
            cocos['sigma_rhothetaphi'] = +1
            cocos['sigma_q_pos'] = +1
            cocos['sign_pprime_pos'] = -1
        elif nbr in [3, 13]:
            cocos['sigma_Bp'] = -1
            cocos['sigma_RphiZ'] = +1
            cocos['sigma_rhothetaphi'] = -1
            cocos['sigma_q_pos'] = -1
            cocos['sign_pprime_pos'] = +1
        elif nbr in [4, 14]:
            cocos['sigma_Bp'] = -1
            cocos['sigma_RphiZ'] = -1
            cocos['sigma_rhothetaphi'] = -1
            cocos['sigma_q_pos'] = -1
            cocos['sign_pprime_pos'] = +1
        elif nbr in [5, 15]:
            cocos['sigma_Bp'] = +1
            cocos['sigma_RphiZ'] = +1
            cocos['sigma_rhothetaphi'] = -1
            cocos['sigma_q_pos'] = -1
            cocos['sign_pprime_pos'] = -1
        elif nbr in [6, 16]:
            cocos['sigma_Bp'] = +1
            cocos['sigma_RphiZ'] = -1
            cocos['sigma_rhothetaphi'] = -1
            cocos['sigma_q_pos'] = -1
            cocos['sign_pprime_pos'] = -1
        elif nbr in [7, 17]:
            cocos['sigma_Bp'] = -1
            cocos['sigma_RphiZ'] = +1
            cocos['sigma_rhothetaphi'] = +1
            cocos['sigma_q_pos'] = +1
            cocos['sign_pprime_pos'] = +1
        elif nbr in [8, 18]:
            cocos['sigma_Bp'] = -1
            cocos['sigma_RphiZ'] = -1
            cocos['sigma_rhothetaphi'] = +1
            cocos['sigma_q_pos'] = +1
            cocos['sign_pprime_pos'] = +1
        else:
            raise ValueError(f"Invalid COCOS number specified: {nbr}.")

        cocos['sign_theta_clockwise'] = cocos['sigma_RphiZ']*cocos['sigma_rhothetaphi']
        return cocos


    def get_flux_surface(self, psi_n, theta=None, closedContourTol = 1e-6):
        """
        Trace the flux surface for the given normalized psi.
        """
        #vertices, _ = self.contour_generator.create_contour(psi_n)
        #v = self.contour_generator.create_contour(psi_n)
        #if len(v) == 1:
        #    vertices = v
        #else:
        #    vertices = v[-1]
            
        vertices = self.contour_generator.create_contour(psi_n)

        if type(vertices) == tuple:
            vertices = vertices[0]

        iClosedContour = None
        for i in range(len(vertices)):
            if len(vertices[i].shape)==2:
                if np.sqrt((vertices[i][0,0]-vertices[i][-1,0])**2 + (vertices[i][0,1]-vertices[i][-1,1])**2)<closedContourTol:
                    if np.amin(vertices[i][:,0]) < self.opoint[0] and np.amax(vertices[i][:,0]) > self.opoint[0] and np.amin(vertices[i][:,1]) < self.opoint[1] and np.amax(vertices[i][:,1]) > self.opoint[1]:
                        iClosedContour = i

        if iClosedContour is not None:
            R, Z = vertices[iClosedContour][:,0], vertices[iClosedContour][:,1]
        else:
            raise ValueError('No closed flux surface was found for psi_n={}'.format(psi_n))

        if theta is not None:
            _theta = np.arctan2(R-self.R0, Z-self.Z0)
            #_theta[-1] = _theta[0] + 2*np.pi
            i = -1
            while _theta[i] < 0:
                _theta[i] += 2*np.pi
                i -= 1

            for i in range(1, _theta.size):
                # The contour finding routine may sometimes give us the
                # same point multiple times, so we have to remove them
                # manually...
                if _theta[i] == _theta[i-1]:
                    _tt = np.zeros((_theta.size-1,))
                    _tt[:i] = _theta[:i]
                    _tt[i:] = _theta[(i+1):]

                    _r = np.zeros((_theta.size-1,))
                    _r[:i] = R[:i]
                    _r[i:] = R[(i+1):]

                    _z = np.zeros((_theta.size-1,))
                    _z[:i] = Z[:i]
                    _z[i:] = Z[(i+1):]

                    _theta = _tt
                    R = _r
                    Z = _z
                    break

            _R = CubicSpline(_theta, R, bc_type='periodic')
            _Z = CubicSpline(_theta, Z, bc_type='periodic')

            R = _R(theta+np.pi/2)
            Z = _Z(theta+np.pi/2)

        return R, Z


    def parametrize_equilibrium(self, psi_n=None, npsi=40):
        """
        Calculates the magnetic equilibrum parameters used by the
        analytical magnetic field in DREAM for a range of flux surfaces.

        :param psi_n: List of normalized poloidal flux for which to calculate the parameters.
        :param npsi:  Number of psi points to calculate parameters for (uniformly spaced between (0, 1]).
        """
        if psi_n is None:
            psi_n = np.linspace(0, 1, npsi+1)[1:]

        radius, psi, kappa, delta, Delta, GOverR0, R, Z = [], [], [], [], [], [], [], []
        for p in psi_n:
            params = self._get_eq_parameters(p)

            radius.append(params['r_minor'])
            psi.append(params['psi'])
            kappa.append(params['kappa'])
            delta.append(params['delta'])
            Delta.append(params['Delta'])
            GOverR0.append(params['GOverR0'])
            R.append(params['R'])
            Z.append(params['Z'])

        radius = np.array(radius)
        psi = np.array(psi) * 2*np.pi / self.R0
        kappa = np.array(kappa)
        delta = np.array(delta)
        Delta = np.array(Delta)
        GOverR0 = np.array(GOverR0)

        return {
            'R0': self.R0,
            'r': radius,
            'psi': psi,
            'kappa': kappa,
            'delta': delta,
            'Delta': Delta,
            'GOverR0': GOverR0
        }


    def _get_eq_parameters(self, psi_n):
        """
        Calculates the magnetic equilibrium parameters used by the
        analytical magnetic field in DREAM for a *SINGLE* flux surface.
        """
        R, Z = self.get_flux_surface(psi_n)

        rho = self.rho(psi_n)
        r_minor = rho * self.a_minor
        Zind = np.argmax(abs(Z))

        R_upper = R[Zind]
        drho_dpsi = self.rho.derivative()(psi_n)
        R0 = self.R0
        R_major = self.R_major(psi_n)

        # Shaping parameters
        psi = self.psi(R[0], Z[0])[0,0]
        kappa = (max(Z)-min(Z)) / (2*r_minor)
        delta = (R_major - R_upper) / r_minor
        #Delta = self.R_major.derivative()(psi_n) / drho_dpsi / self.a_minor
        Delta = (max(R)+min(R))/2 - R0
        GOverR0 = self.f_psi(psi_n) / R0

        return {
            'r_minor': r_minor,
            'psi': psi,
            'kappa': kappa,
            'delta': delta,
            'Delta': Delta,
            'GOverR0': GOverR0,
            'R': R,
            'Z': Z
        }


    def get_Bpol_coeff(self):
        """
        Return the COCOS conversion coefficient for the poloidal
        magnetic field.
        """
        cocos = self.cocos()
        coeff = cocos['sigma_RphiZ']*cocos['sigma_Bp']/(2*np.pi)**cocos['exp_Bp']
        return coeff
        
    def get_Br(self, R, Z):
        """
        Return the radial magnetic field component on the given (R,Z) grid.
        """
        coeff = self.get_Bpol_coeff()
        Br = coeff/R * self.psi(R, Z, dy=1, grid=False)
        return Br


    def get_Bz(self, R, Z):
        """
        Return the vertical magnetic field component on the given (R,Z) grid.
        """
        coeff = self.get_Bpol_coeff()
        Bz = -coeff/R * self.psi(R, Z, dx=1, grid=False)
        return Bz


    def get_Btor(self, R, Z):
        """
        Return the toroidal magnetic field component on the given (R,Z) grid.
        """
        psi = self.psi(R, Z, grid=False)
        psi_n = (psi - self.psi_axis) / (self.psi_bdry - self.psi_axis)

        f  = self.f_psi(psi_n)
        Btor = f / R

        return Btor


    def getBOfTheta(self, psi_n):
        """
        Evaluates the magnetic field strength along the specified flux
        surface as a function of poloidal angle.
        """
        R, Z = self.get_flux_surface(psi_n)

        Br = self.get_Br(R, Z)
        Bz = self.get_Bz(R, Z)
        Bt = self.get_Btor(R, Z)

        B = np.sqrt(Br**2 + Bz**2 + Bt**2)
        theta = np.arctan2(R-self.R0, Z-self.Z0)

        return theta, B


    def load(self, filename, override_psilim=False):
        """
        Load data from the named GEQDSK file to this object.
        """
        self.process_data(data, override_psilim=override_psilim)


    def process_data(self, data, override_psilim=False):
        """
        Load data from the given EQDSK dictionary.
        """
        self.nr = data['nx']
        self.nz = data['ny']
        self.cocos_number = data['cocos']

        psi = data['psi']
        self.bcentr = data['bcentr']
        self.psi_axis = data['psiaxis']
        self.psi_bdry = data['psiedge']

        psi_n = np.linspace(0, 1, self.nr)

        self.f_psi    = InterpolatedUnivariateSpline(psi_n, data["fpol"])
        self.ff_prime = InterpolatedUnivariateSpline(psi_n, data["ffprime"])
        self.q        = InterpolatedUnivariateSpline(psi_n, data["q"])
        self.pressure = InterpolatedUnivariateSpline(psi_n, data["p"])
        self.p_prime  = self.pressure.derivative()

        self.R0 = data['raxis']
        self.Z0 = data['zaxis']
        self.R = np.linspace(data["rleft"], data["rleft"]+data["rboxlen"], self.nr)
        self.Z = np.linspace(data["zmid"]-data["zboxlen"]/2, data["zmid"]+data["zboxlen"]/2, self.nz)

        self.psi = RectBivariateSpline(self.R, self.Z, psi)

        if 'rlim' in  data and 'zlim' in data:
            self.rlim = data['rlim']
            self.zlim = data['zlim']
        else:
            self.rlim = []
            self.zlim = []

        if 'rplas' in data and 'zplas' in data:
            self.rplas = data['rplas']
            self.zplas = data['zplas']
        else:
            self.rplas = []
            self.zplas = []

        # Locate O-point
        self.opoint = self.find_o_point()

        if override_psilim:
            self.psi_axis = self.psi(*self.opoint)
            # The -1e-4 term is required for us to be able to find all psi_n on [0, 1]
            self.psi_bdry = self.psi(self.rplas[0], self.zplas[0]) - 2e-4

        # Set up contour generator
        psi2d = np.transpose(self.psi(self.R, self.Z))
        psin2d = (psi2d - self.psi_axis) / (self.psi_bdry - self.psi_axis)
        R, Z = np.meshgrid(self.R, self.Z)
        self.contour_generator = QuadContourGenerator(R, Z, psin2d, None, True, 0)

        rho = np.zeros(psi_n.shape)
        R_major = np.zeros(psi_n.shape)

        for i, i_psiN in enumerate(psi_n[1:]):
            surface_R, surface_Z = self.get_flux_surface(psi_n=i_psiN)

            rho[i+1] = (max(surface_R)-min(surface_R)) / 2
            R_major[i+1] = (max(surface_R)+min(surface_R)) / 2

        self.lcfs_R = surface_R
        self.lcfs_Z = surface_Z

        self.a_minor = rho[-1]

        rho = rho / rho[-1]

        R_major[0] = R_major[1] + psi_n[1] * (R_major[2]-R_major[1]) / (psi_n[2]-psi_n[1])

        self.rho = InterpolatedUnivariateSpline(psi_n, rho)
        self.R_major = InterpolatedUnivariateSpline(psi_n, R_major)

        self.R0 = self.R_major(0)


    def plot_flux_surfaces(self, ax=None, nr=10, ntheta=200, fit=True, *args, **kwargs):
        """
        Plot the flux surfaces of this magnetic equilibrium.

        :param ax:       Matplotlib Axes object to use for drawing.
        :param nr:       Number of flux surfaces to plot.
        :param ntheta:   Number of poloidal angles to plot (for contour fit).
        :param fit:      If ``True``, plots the DREAM parameter fit surfaces instead of the actual flux surfaces.
        :param *args:    Arguments for ``ax.plot()``.
        :param **kwargs: Keyword arguments for ``ax.plot()``.
        """
        fig = None
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        if fit:
            theta = np.linspace(0, 2*np.pi, ntheta)
            p = self.parametrize_equilibrium(npsi=nr)

            for i in range(p['r'].size):
                R = p['R0'] + p['Delta'][i] + p['r'][i]*np.cos(theta + p['delta'][i]*np.sin(theta))
                Z = p['r'][i]*p['kappa'][i]*np.sin(theta)

                ax.plot(R, Z, *args, **kwargs)
            ax.axis('equal')
        else:
            psi_n = np.linspace(0, 1, nr+1)[1:]
            for p in psi_n:
                R, Z = self.get_flux_surface(p)
                ax.plot(R, Z, *args, **kwargs)
            ax.axis('equal')

        if len(self.rlim) > 0:
            ax.plot(self.rlim, self.zlim, 'k', linewidth=2)
        if len(self.rplas) > 0:
            ax.plot(self.rplas, self.zplas, 'b', linewidth=2)

        ax.plot(self.R0, self.Z0, 'rd', linewidth=2)

        return ax


    def get_LUKE(self, npsi=80, ntheta=90):
        """
        Returns equilibrium data in the LUKE equilibrium format.
        """
        theta = np.linspace(0, 2*np.pi, ntheta)
        psi_n = np.linspace(0, 1, npsi+1)[1:]

        Rp, Zp = self.R0, self.Z0
        psi_apRp = 2*np.pi * self.psi(Rp+self.rho(psi_n), self.Z0) * self.a_minor / Rp

        ptx = np.zeros((psi_n.size, ntheta))
        pty = np.zeros((psi_n.size, ntheta))
        for i in range(npsi):
            ptx[i,:], pty[i,:] = self.get_flux_surface(psi_n[i], theta=theta)


        ptBx = self.get_Br(ptx, pty)
        ptBy = self.get_Bz(ptx, pty)
        ptBPHI = self.get_Btor(ptx, pty)

        return {
            'id': 'GEQDSK data',
            'Rp': np.array([Rp]), 'Zp': np.array([Zp]),
            'psi_apRp': psi_apRp,
            'theta': theta,
            'ptx': ptx.T-Rp, 'pty': pty.T-Zp,
            'ptBx': ptBx.T, 'ptBy': ptBy.T, 'ptBPHI': ptBPHI.T
        }


    def find_o_point(self):
        """
        Locate the (R, Z) coordinates of the O-point.
        """
        cocos = self.cocos()
        sBp = cocos['sigma_Bp']

        if scipy.__version__ >= '1.9.0':
            def jac(x):
                dpsi_dx = self.psi.partial_derivative(1, 0)
                dpsi_dy = self.psi.partial_derivative(0, 1)

                return np.array([dpsi_dx(*x), dpsi_dy(*x)]).flatten()
        else:
            jac = None

        res = scipy.optimize.minimize(lambda x : self.psi(*x)*sBp, x0=(self.R0, self.Z0), jac=jac)

        if res.success:
            return res.x
        else:
            raise Exception(f"Unable to locate O-point: {res.message}.")


    def get_SOFT(self, nr=80, nz=90, find_maxis=True):
        """
        Returns equilibrium data in the SOFT equilibrium format.
        """
        if len(self.rplas) == 0 or len(self.rlim) == 0:
            raise Exception("Generation of SOFT equilibrium requires a wall and plasma boundary to be specified.")

        rmin, rmax = np.amin(self.rlim), np.amax(self.rlim)
        zmin, zmax = np.amin(self.zlim), np.amax(self.zlim)
        #R, Z = np.meshgrid(np.linspace(rmin, rmax, nr), np.linspace(zmin, zmax, nz))
        r, z = np.linspace(rmin, rmax, nr), np.linspace(zmin, zmax, nz)
        R, Z = np.meshgrid(r, z)

        verBr = self.get_Br(R[0,:], Z[0,:])
        verBz = self.get_Bz(R[0,:], Z[0,:])
        verBphi = self.get_Btor(R[0,:], Z[0,:])

        raxis, zaxis = self.opoint

        return {
            'name': 'EQDSK',
            'desc': 'EQDSK',
            'cocos': self.cocos_number,
            'maxis': np.array([raxis, zaxis]),
            'Br': self.get_Br(R, Z),
            'Bz': self.get_Bz(R, Z),
            'Bphi': self.get_Btor(R, Z),
            'Psi': self.psi(r, z),
            'verBr': verBr,
            'verBz': verBz,
            'verBphi': verBphi,
            'separatrix': np.array([self.rplas, self.zplas]),
            'wall': np.array([self.rlim, self.zlim]),
            'r': r, 'z': z,
        }


    def save_eq_parameters(self, filename, nr=40):
        """
        Save the DREAM analytical equilibrium parameters corresponding to this
        GEQDSK file to an HDF5 file named ``filename``.
        """
        params = self.parametrize_equilibrium(npsi=nr)

        with h5py.File(filename, 'w') as f:
            f['r'] = params['r']
            f['Delta'] = params['Delta']
            f['delta'] = params['delta']
            f['GOverR0'] = params['GOverR0']
            f['kappa'] = params['kappa']
            f['psi_p'] = params['psi']

    
    def save_LUKE(self, filename, npsi=80, ntheta=90):
        """
        Save this equilibrium in a LUKE compatible equilibrium file.
        """
        equil = self.get_LUKE(npsi=npsi, ntheta=ntheta)

        with h5py.File(filename, 'w') as f:
            f.create_group('equil')

            for key in equil.keys():
                f[f'equil/{key}'] = equil[key]


    def save_SOFT(self, filename, nr=80, nz=90, find_maxis=True):
        """
        Save this equilibrium in a SOFT compatible equilibrium file.
        """
        equil = self.get_SOFT(nr=nr, nz=nz, find_maxis=find_maxis)

        with h5py.File(filename, 'w') as f:
            for key in equil.keys():
                f[key] = equil[key]


