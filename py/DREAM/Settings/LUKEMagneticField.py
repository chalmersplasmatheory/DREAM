# Implementation of LUKE magnetic equilibrium data

import h5py
import matplotlib.pyplot as plt
import numpy as np
from .. import DREAMIO
from .. import helpers
from . NumericalMagneticField import NumericalMagneticField


class LUKEMagneticField(NumericalMagneticField):
    

    def __init__(self, filename):
        """
        Constructor.
        """
        self.load(filename)

        super().__init__(a=self.getMinorRadius(), R0=self.getMajorRadius())


    def load(self, filename, path=''):
        """
        Load a LUKE magnetic equilibrium from the named file.
        """
        PATH = '{}/equil'.format(path)

        with h5py.File(filename, 'r') as f:
            self.id = DREAMIO.getData(f[PATH], 'id')
            if type(self.id) != str:
                self.id = ''

            self.Rp = f['{}/Rp'.format(PATH)][:]
            self.Zp = f['{}/Zp'.format(PATH)][:]
            self.psi_apRp = f['{}/psi_apRp'.format(PATH)][:]
            self.theta = f['{}/theta'.format(PATH)][:]
            self.ptx = f['{}/ptx'.format(PATH)][:]
            self.pty = f['{}/pty'.format(PATH)][:]
            self.ptBx = f['{}/ptBx'.format(PATH)][:]
            self.ptBy = f['{}/ptBy'.format(PATH)][:]
            self.ptBPHI = f['{}/ptBPHI'.format(PATH)][:]

        if self.Rp.ndim == 2:
            self.Rp = self.Rp[0,0]
        if self.Zp.ndim == 2:
            self.Zp = self.Zp[0,0]
        if self.psi_apRp.shape[0] == 1:
            self.psi_apRp = self.psi_apRp[0,:]
        if self.theta.shape[0] == 1:
            self.theta = self.theta[0,:]


    def getMinorRadius(self):
        """
        Returns plasma minor radius.
        """
        return float(self.ptx[0,-1])


    def getMajorRadius(self):
        """
        Returns tokamak major radius.
        """
        return float(self.Rp)


    def visualize(self, npsi=20, ax=None, show=None):
        """
        Visualize this magnetic field.
        """
        red   = (249/255, 65/255, 68/255)
        black = (87/255, 117/255, 144/255)
        gray  = (120/255, 120/255, 120/255)

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True
    
        R, Z = self.ptx[:,:-1]+self.Rp, self.pty[:,:-1]+self.Zp
        if npsi > R.shape[1]:
            npsi = R.shape[1]

        sel  = [int(i) for i in np.round(np.linspace(0, R.shape[1]-1, npsi))]

        ax.plot(R[:,sel], Z[:,sel], color=gray, linewidth=0.5)
        ax.plot(self.ptx[:,-1] + self.Rp, self.pty[:,-1] + self.Zp, color=red, linewidth=2)
        ax.plot(self.Rp, self.Zp, 's', color=red)

        ax.set_title(helpers.safeTeXstring(self.id))
        ax.set_xlabel('$R$ (m)')
        ax.set_ylabel('$Z$ (m)')

        ax.axis('equal')

        if show:
            plt.show()

        return ax

