# Implementation of LUKE magnetic equilibrium data

import h5py
import matplotlib.pyplot as plt
import numpy as np
from desc.grid import LinearGrid
import desc.io
from .. import DREAMIO
from .. import helpers
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
        self.B = None
        self.BdotGradPhi = None
        self.Jacobian = None
        self.g_tt = None
        self.g_tp = None
        self.lambda_t = None
        self.lambda_p = None

        self.eq = desc.io.load(self.filename)

        self.grid = LinearGrid(L=self.nr - 1, M=int((self.ntheta - 1) / 2), N=int((self.nphi - 1) / 2), endpoint=True, NFP=self.eq.NFP)

        self.R0 = self.eq.compute('R0', grid=self.grid)['R0']
        self.a = self.eq.compute('a', grid=self.grid)['a']

        self.rho = self.grid.nodes[self.grid.unique_rho_idx,0]
        self.theta = self.grid.nodes[self.grid.unique_theta_idx,1]
        self.phi = self.grid.nodes[self.grid.unique_zeta_idx,2]



    def load(self):
        """
        Load a DESC magnetic equilibrium from the named file.
        """

        self.f_passing = 1- self.eq.compute('trapped fraction', grid=self.grid)['trapped fraction'][self.grid.unique_rho_idx]
        self.B_min = self.eq.compute('min_tz |B|', grid=self.grid)['min_tz |B|'][self.grid.unique_rho_idx]
        self.B_max = self.eq.compute('max_tz |B|', grid=self.grid)['max_tz |B|'][self.grid.unique_rho_idx]

        self.G = self.eq.compute('G', grid=self.grid)['G'][self.grid.unique_rho_idx]
        self.I = self.eq.compute('I', grid=self.grid)['I'][self.grid.unique_rho_idx]
        iota = self.eq.compute('iota', grid=self.grid)['iota']
        self.iota = iota[self.grid.unique_rho_idx]

        self.B = self.eq.compute('|B|', grid=self.grid)['|B|']
        self.BdotGradPhi = self.eq.compute('B^zeta', grid=self.grid)['B^zeta']
        self.Jacobian = self.eq.compute('sqrt(g)', grid=self.grid)['sqrt(g)']

        self.g_tt = self.eq.compute('g_tt', grid=self.grid)['g_tt']
        self.g_tp = self.eq.compute('g_tz', grid=self.grid)['g_tz']
        self.lambda_t = self.eq.compute('lambda_t', grid=self.grid)['lambda_t'] + iota*self.eq.compute('nu_t', grid=self.grid)['nu_t']
        self.lambda_p = self.eq.compute('lambda_z', grid=self.grid)['lambda_z'] + iota*self.eq.compute('nu_z', grid=self.grid)['nu_z']



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


    def visualize(self, phi, nrho=None, ntheta=None, ax=None, show=None): # TODO
        """
        Visualize this magnetic field.
        """
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
            nrho = self.nr
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

