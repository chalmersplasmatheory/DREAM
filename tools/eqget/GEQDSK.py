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


class GEQDSK:
    

    def __init__(self, filename):
        """
        Constructor.
        """
        self.load(filename)


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


    def get_flux_surface(self, psi_n, theta=None):
        """
        Trace the flux surface for the given normalized psi.
        """
        vertices, _ = self.contour_generator.create_contour(psi_n)
        
        R, Z = vertices[0][:,0], vertices[0][:,1]

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


    def get_Br(self, R, Z):
        """
        Return the radial magnetic field component on the given (R,Z) grid.
        """
        Br = 1/R * self.psi(R, Z, dy=1, grid=False)
        return Br


    def get_Bz(self, R, Z):
        """
        Return the vertical magnetic field component on the given (R,Z) grid.
        """
        Bz = 1/R * self.psi(R, Z, dx=1, grid=False)
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


    def load(self, filename):
        """
        Load data from the named GEQDSK file to this object.
        """
        data = self.load_geqdsk(filename)
        self.process_data(data)


    def load_geqdsk(self, filename, cocos=1):
        """
        Load the named GEQDSK file.
        """
        with open(filename) as fh:
            header = fh.readline()
            words = header.split()
            if len(words) < 3:
                raise ValueError("Expected at least 3 numbers on first line")

            nx, ny = int(words[-2]), int(words[-1])
            
            data = {"nx": nx, "ny": ny}
            fields = ["rdim", "zdim", "rcentr", "rleft", "zmid", "rmagx",
                      "zmagx", "simagx", "sibdry", "bcentr", "cpasma", "simagx",
                      None, "rmagx", None, "zmagx", None, "sibdry", None, None]

            values = self._next_value(fh)
            
            for f in fields:
                val = next(values)
                if f:
                    data[f] = val

            def _read_1d(n):
                """
                Read a 1D array of length n from the GEQDSK file.
                """
                val = np.zeros(n)
                for i in range(n):
                    val[i] = next(values)

                return val


            def _read_2d(n, m):
                """
                Read a 2D (n,m) array in Fortran order
                """
                val = np.zeros((n, m))
                for j in range(m):
                    for i in range(n):
                        val[i, j] = next(values)

                return val


            data["fpol"] = _read_1d(nx)
            data["pres"] = _read_1d(nx)
            data["ffprime"] = _read_1d(nx)
            data["pprime"] = _read_1d(nx)

            data["psi"] = _read_2d(nx, ny)

            data["qpsi"] = _read_1d(nx)

            # Ensure that psi is divided by 2pi
            if cocos > 10:
                for var in ["psi", "simagx", "sibdry"]:
                    data[var] /= 2 * pi

            nbdry = next(values)
            nlim = next(values)

            if nbdry > 0:
                data["rbdry"] = np.zeros(nbdry)
                data["zbdry"] = np.zeros(nbdry)
                for i in range(nbdry):
                    data["rbdry"][i] = next(values)
                    data["zbdry"][i] = next(values)

            if nlim > 0:
                data["rlim"] = np.zeros(nlim)
                data["zlim"] = np.zeros(nlim)
                for i in range(nlim):
                    data["rlim"][i] = next(values)
                    data["zlim"][i] = next(values)

            return data


    def process_data(self, data):
        """
        Load data from the given GEQDSK dictionary.
        """
        self.nr = data['nx']
        self.nz = data['ny']

        psi = data['psi']
        self.bcentr = data['bcentr']
        self.psi_axis = data['simagx']
        self.psi_bdry = data['sibdry']

        psi_n = np.linspace(0, 1, self.nr)

        self.f_psi    = InterpolatedUnivariateSpline(psi_n, data["fpol"])
        self.ff_prime = InterpolatedUnivariateSpline(psi_n, data["ffprime"])
        self.q        = InterpolatedUnivariateSpline(psi_n, data["qpsi"])
        self.pressure = InterpolatedUnivariateSpline(psi_n, data["pres"])
        self.p_prime  = self.pressure.derivative()

        self.Z0 = data['zmid']
        self.R = np.linspace(data["rleft"], data["rleft"]+data["rdim"], self.nr)
        self.Z = np.linspace(data["zmid"]-data["zdim"]/2, data["zmid"]+data["zdim"]/2, self.nz)

        self.psi = RectBivariateSpline(self.R, self.Z, psi)

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


