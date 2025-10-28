# Implementation of LUKE magnetic equilibrium data

import h5py
import matplotlib.pyplot as plt
import numpy as np
#from .. import DREAMIO
#from .. import helpers
from . NumericalMagneticField import NumericalMagneticField




class StellaratorMagneticField(NumericalMagneticField):
    

    def __init__(self, filename, nr, ntheta=129, nphi=129): # TODO: OK?
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

        self.grid = LinearGrid(L=self.nr - 1, M=int((self.ntheta - 1) / 2), N=int((self.nphi - 1) / 2), endpoint=True, NFP=self.eq.NFP)

        R = self.eq.compute('R', grid=self.grid)['R'].reshape((self.nphi, self.nr, self.ntheta))
        self.R0 = float(R[0,0,0])
        self.a = float(R[0,-1,0]-R[0,0,0])
        #self.R0 = float(self.eq.compute('R0', grid=self.grid)['R0']) # TODO: This is not the major radius as defined in DREAM
        #self.a = float(self.eq.compute('a', grid=self.grid)['a']) # TODO: This is not the minor radius as defined in DREAM
        self.nfp = int(self.eq.NFP)

        self.rho = np.array(self.grid.nodes[self.grid.unique_rho_idx,0], dtype=np.float64)
        self.theta = np.array(self.grid.nodes[self.grid.unique_theta_idx,1], dtype=np.float64)
        self.phi = np.array(self.grid.nodes[self.grid.unique_zeta_idx,2], dtype=np.float64)


    def load(self, savedata=False): # TODO: Spara behandlad data till fil
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
        self.Jacobian = np.array(self.eq.compute('sqrt(g)', grid=self.grid)['sqrt(g)'], dtype=np.float64)

        self.g_tt = np.array(self.eq.compute('g_tt', grid=self.grid)['g_tt'], dtype=np.float64)
        self.g_tp = np.array(self.eq.compute('g_tz', grid=self.grid)['g_tz'], dtype=np.float64)
        self.lambda_t = np.array(self.eq.compute('lambda_t', grid=self.grid)['lambda_t'] + iota*self.eq.compute('nu_t', grid=self.grid)['nu_t'], dtype=np.float64)
        self.lambda_p = np.array(self.eq.compute('lambda_z', grid=self.grid)['lambda_z'] + iota*self.eq.compute('nu_z', grid=self.grid)['nu_z'], dtype=np.float64)



        if savedata:
            dic = { "R0" : self.R0,
                    "a": self.a,
                    "nr": self.nr,
                    "nfp": self.nfp,
                    "rho": self.theta,
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
            with h5py.File(f"{self.savefilename}.h5", 'w') as hf:
                for key, value in dic.items():
                    hf[key] = value


        from desc.grid import LinearGrid
        # TODO: Temporary for benchmark, remove
        lukegrid = LinearGrid(L=self.nr-1, M=int((self.ntheta - 1) / 2), zeta=0, endpoint=True)
        theta = lukegrid.nodes[lukegrid.unique_theta_idx, 1]
        R = self.eq.compute('R', grid=lukegrid)['R'].reshape((self.nr, self.ntheta)).T
        R -= R[0,0]
        psi = 2 * np.pi * self.eq.compute('chi', grid=lukegrid)['chi'][lukegrid.unique_rho_idx] / self.R0 * self.a

        Z = self.eq.compute('Z', grid=lukegrid)['Z'].reshape((self.nr, self.ntheta)).T
        B_R = self.eq.compute('B_R', grid=lukegrid)['B_R'].reshape((self.nr, self.ntheta)).T
        B_Z = self.eq.compute('B_Z', grid=lukegrid)['B_Z'].reshape((self.nr, self.ntheta)).T
        B_phi = self.eq.compute('B_phi', grid=lukegrid)['B_phi'].reshape((self.nr, self.ntheta)).T

        dictemp = {#"/equil":
                        #{
                         "/equil/id": 1,
                         "/equil/Rp": np.array([self.R0]),
                         "/equil/Zp": np.array([0]),
                         "/equil/psi_apRp": psi.reshape((-1,1)),
                         "/equil/theta": theta,
                         "/equil/ptx": R,
                         "/equil/pty": Z,
                         "/equil/ptBx": B_R,
                         "/equil/ptBy": B_Z,
                         "/equil/ptBPHI": B_phi
                        #}
        }

        with h5py.File(f"data/numericalLUKE.h5", 'w') as hf:
            for key, value in dictemp.items():
                hf[key] = value
        '''
        with h5py.File(f"numericalLUKE.h5", 'w') as hf:
            for group_name, group_data in dictemp.items():
                group = hf.create_group(group_name)
                for key, value in group_data.items():
                    if isinstance(value, np.ndarray):
                        group.create_dataset(key, data=value)
                    else:
                        group.attrs[key] = value
        '''
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
            #ax.plot(R[0,0], Z[0,0], 's', color=red)

        ax.plot(self.R0, 0, 's', color=red)#'X', color=black)

        ax.set_xlabel('$R$ (m)')
        ax.set_ylabel('$Z$ (m)')

        ax.axis('equal')

        if show:
            plt.show()

        return ax

