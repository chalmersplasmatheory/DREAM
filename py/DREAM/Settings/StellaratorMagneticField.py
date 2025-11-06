# Implementation of LUKE magnetic equilibrium data

import h5py
import matplotlib.pyplot as plt
import numpy as np
#from .. import DREAMIO
#from .. import helpers
from . NumericalMagneticField import NumericalMagneticField




class StellaratorMagneticField(NumericalMagneticField):
    

    def __init__(self, filename, nr, ntheta=129, nphi=129, datafilename=None): # TODO: OK?
        """
        Constructor.
        """
        self.filename = filename
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi

        self.rho = None
        self.theta = None
        self.phi = None
        self.f_passing = None
        self.B_min = None
        self.B_max = None
        self.G = None
        self.I = None
        self.iota = None
        self.psi_T = None
        self.B = None
        self.BdotGradPhi = None
        self.Jacobian = None
        self.g_tt = None
        self.g_tp = None
        self.lambda_t = None
        self.lambda_p = None

        from desc.grid import LinearGrid
        import desc.io

        self.eq = desc.io.load(self.filename)

        if datafilename is None:
            self.grid = LinearGrid(L=int(self.nr-1), M=int((self.ntheta - 1) / 2), N=int((self.nphi - 1) / 2), endpoint=True, NFP=self.eq.NFP)

            self.R = np.array(self.eq.compute('R', grid=self.grid)['R'], dtype=np.float64)
            self.Z = np.array(self.eq.compute('Z', grid=self.grid)['Z'], dtype=np.float64)
            self.R0 = float(self.eq.axis.R_n[self.eq.axis.R_basis.get_idx(0)])
            self.a = float(self.eq.compute('a', grid=self.grid)['a'])
            self.nfp = int(self.eq.NFP)

            self.rho = np.array(self.grid.nodes[self.grid.unique_rho_idx,0]*self.a, dtype=np.float64)
            self.theta = np.array(self.grid.nodes[self.grid.unique_theta_idx,1], dtype=np.float64)
            self.phi = np.array(self.grid.nodes[self.grid.unique_zeta_idx,2], dtype=np.float64)
        else:
            with h5py.File(f"{datafilename}.h5", 'r') as hf:
                self.R = np.array(hf['R'][:], dtype=np.float64)
                self.Z = np.array(hf['Z'][:], dtype=np.float64)
                self.R0 = hf["R0"][()]
                self.a = hf["a"][()]
                self.nfp = hf["nfp"][()]
                self.rho = np.array(hf["rho"][:], dtype=np.float64)
                self.theta = np.array(hf["theta"][:], dtype=np.float64)
                self.phi =np.array( hf["phi"][:], dtype=np.float64)
                self.f_passing = np.array(hf["f_passing"][:], dtype=np.float64)
                self.B_min = np.array(hf["B_min"][:], dtype=np.float64)
                self.B_max = np.array(hf["B_max"][:], dtype=np.float64)
                self.G = np.array(hf["G"][:], dtype=np.float64)
                self.I = np.array(hf["I"][:], dtype=np.float64)
                self.iota = np.array(hf["iota"][:], dtype=np.float64)
                self.psi_T = np.array(hf["psi_T"][:], dtype=np.float64)
                self.B = np.array(hf["B"][:], dtype=np.float64)
                self.BdotGradPhi = np.array(hf["BdotGradPhi"][:], dtype=np.float64)
                self.Jacobian = np.array(hf["Jacobian"][:], dtype=np.float64)
                self.g_tt = np.array(hf["g_tt"][:], dtype=np.float64)
                self.g_tp = np.array(hf["g_tp"][:], dtype=np.float64)
                self.lambda_t = np.array(hf["lambda_t"][:], dtype=np.float64)
                self.lambda_p = np.array(hf["lambda_p"][:], dtype=np.float64)


    def load(self, savedata=True, savefilename='numericStellaratorSettings'): # TODO: Spara behandlad data till fil
        """
        Load a DESC magnetic equilibrium from the named file.
        """

        self.f_passing = np.array(1- self.eq.compute('trapped fraction', grid=self.grid)['trapped fraction'][self.grid.unique_rho_idx], dtype=np.float64)
        self.B_min = np.array(self.eq.compute('min_tz |B|', grid=self.grid)['min_tz |B|'][self.grid.unique_rho_idx], dtype=np.float64)
        self.B_max = np.array(self.eq.compute('max_tz |B|', grid=self.grid)['max_tz |B|'][self.grid.unique_rho_idx], dtype=np.float64)
        self.G = np.array(self.eq.compute('G', grid=self.grid)['G'][self.grid.unique_rho_idx], dtype=np.float64)
        self.I = np.array(self.eq.compute('I', grid=self.grid)['I'][self.grid.unique_rho_idx], dtype=np.float64)
        iota = np.array(self.eq.compute('iota', grid=self.grid)['iota'], dtype=np.float64)
        self.iota = iota[self.grid.unique_rho_idx]
        self.psi_T = np.array(self.eq.compute('Psi', grid=self.grid)['Psi'][self.grid.unique_rho_idx], dtype=np.float64)

        self.B = np.array(self.eq.compute('|B|', grid=self.grid)['|B|'], dtype=np.float64)
        self.BdotGradPhi = np.array(self.eq.compute('B^zeta', grid=self.grid)['B^zeta'], dtype=np.float64)
        self.Jacobian = np.array(self.eq.compute('sqrt(g)', grid=self.grid)['sqrt(g)'], dtype=np.float64) / self.a
        self.g_tt = np.array(self.eq.compute('g_tt', grid=self.grid)['g_tt'], dtype=np.float64)
        self.g_tp = np.array(self.eq.compute('g_tz', grid=self.grid)['g_tz'], dtype=np.float64)
        isAxis = np.where(self.grid.nodes[:,0] == 0)
        self.g_tt = self.g_tt / self.Jacobian**2
        self.g_tp = self.g_tp / self.Jacobian**2
        for i in isAxis:
            self.g_tt[i] = self.g_tt[i+self.ntheta]
            self.g_tp[i] = self.g_tp[i+self.ntheta]

        self.lambda_t = np.array(self.eq.compute('lambda_t', grid=self.grid)['lambda_t'] + iota*self.eq.compute('nu_t', grid=self.grid)['nu_t'], dtype=np.float64)
        self.lambda_p = np.array(self.eq.compute('lambda_z', grid=self.grid)['lambda_z'] + iota*self.eq.compute('nu_z', grid=self.grid)['nu_z'], dtype=np.float64)

        if savedata:
            dic = { "R" : self.R,
                    "Z" : self.Z,
                    "R0" : self.R0,
                    "a": self.a,
                    "nfp": self.nfp,
                    "rho": self.rho,
                    "theta": self.theta,
                    "phi": self.phi,
                    "f_passing": self.f_passing,
                    "B_min": self.B_min,
                    "B_max": self.B_max,
                    "G": self.G,
                    "I": self.I,
                    "iota": self.iota,
                    "psi_T": self.psi_T,
                    "B": self.B,
                    "BdotGradPhi": self.BdotGradPhi,
                    "Jacobian": self.Jacobian,
                    "g_tt": self.g_tt,
                    "g_tp": self.g_tp,
                    "lambda_t": self.lambda_t,
                    "lambda_p": self.lambda_p
            }
            with h5py.File(f"{savefilename}.h5", 'w') as hf:
                for key, value in dic.items():
                    hf[key] = value


        with h5py.File(f"data/numericalLUKE.h5", 'w') as hf:
            for key, value in dictemp.items():
                hf[key] = value

        # TODO: Save this as quantities instead of R and Z, for ouput and for SPI in the future
        #R_mn = self.eq.compute('lambda_z', grid=self.grid)['lambda_z'] + iota * \
        #           self.eq.compute('nu_z', grid=self.grid)['nu_z']



    def getMinorRadius(self):
        """
        Returns plasma minor radius.
        """
        return float(self.a)


    def getMajorRadius(self):
        """
        Returns tokamak major radius.
        """
        return float(self.R0)


    def visualize(self, phi, nrho=None, ntheta=None, ax=None, show=None):
        """
        Visualize this magnetic field.
        """

        from desc.grid import LinearGrid
        import desc.io

        red = (249 / 255, 65 / 255, 68 / 255)
        black = (87 / 255, 117 / 255, 144 / 255)
        gray = (120 / 255, 120 / 255, 120 / 255)

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        if np.isscalar(phi):
            phi = [phi]

        if nrho == None:
            nrho = self.nr-1
        if ntheta == None:
            ntheta = self.ntheta

        for iphi in phi:
            plotgrid = LinearGrid(L=nrho, M=int((ntheta - 1) / 2), zeta=iphi, endpoint=True)

            R, Z = self.eq.compute('R', grid=plotgrid)['R'].reshape((nrho+1, ntheta)).T, self.eq.compute('Z', grid=plotgrid)['Z'].reshape((nrho+1, ntheta)).T

            ax.plot(R[:,:-1], Z[:,:-1], color=gray, linewidth=0.5)
            ax.plot(R[:,-1], Z[:,-1], color=red, linewidth=2)
            ax.plot(R[0,0], Z[0,0], 's', color=red)

        ax.plot(self.R0, 0, 'X', color=black)
        ax.set_xlabel('$R$ (m)')
        ax.set_ylabel('$Z$ (m)')

        ax.axis('equal')

        if show:
            plt.show()

        return ax

