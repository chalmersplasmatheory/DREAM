
import matplotlib.pyplot as plt
import numpy as np


class Equilibrium:
    

    def __init__(self, eq=None):
        """
        Constructor.

        :param eq: Equilibrium data from DREAM output.
        """
        if eq is not None:
            self.setEquilibrium(eq)


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        s  = f'TOKAMAK EQUILIBRIUM (ntheta = {self.RMinusR0.shape[0]}, npsi = {self.RMinusR0.shape[1]})\n'
        s += f'  Axis at ({self.R0[0]}, {self.Z0[0]})\n'
        s += f'  Minor radius a = {self.RMinusR0_f[0,-1]} m\n'

        return s


    def setEquilibrium(self, eq):
        """
        Set equilibrium data based on output from DREAM.
        """
        self.R0 = eq['R0']
        self.Z0 = eq['Z0']
        self.RMinusR0 = eq['RMinusR0']
        self.RMinusR0_f = eq['RMinusR0_f']
        self.ZMinusZ0 = eq['ZMinusZ0']
        self.ZMinusZ0_f = eq['ZMinusZ0_f']
        self.theta = eq['theta']


    def visualize(self, shifted=False, maxis=True, ax=None, show=None, **kwargs):
        """
        Visualize this equilibrium.

        :param shifted: If ``True``, shifts the surfaces so that (0, 0) is on the magnetic axis.
        """
        red  = (249/255, 65/255, 68/255)
        black = (87/255, 117/255, 144/255)
        gray = (190/255, 190/255, 190/255)

        # Set up axes (if not already done)

        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        rn, zn = 0, 0
        if not shifted and not np.isinf(self.R0):
            rn = self.R0
            zn = self.Z0

        # Flux surfaces
        ax.plot(self.RMinusR0+rn, self.ZMinusZ0+zn, color=gray, linewidth=1, **kwargs)
        # ...close the flux surfaces
        ax.plot(np.array([self.RMinusR0[-1,:], self.RMinusR0[0,:]])+rn, np.array([self.ZMinusZ0[-1,:], self.ZMinusZ0[0,:]])+zn, color=gray, linewidth=1, **kwargs)
        # Limiter
        ax.plot(self.RMinusR0_f[:,-1]+rn, self.ZMinusZ0_f[:,-1]+zn, color=black, linewidth=2, **kwargs)
        ax.plot(self.RMinusR0_f[(0,-1),-1]+rn, self.ZMinusZ0_f[(0,-1),-1]+zn, color=black, linewidth=2, **kwargs)
        # Magnetic axis 
        if maxis:
            ax.plot(rn, zn, 'x', color=red)
        ax.axis('equal')

        if np.isinf(self.R0):
            ax.set_xlabel('Major radius $R-R_0$ (m)')
        else:
            ax.set_xlabel('Major radius $R$ (m)')
        ax.set_ylabel('Height $Z$ (m)')

        if show:
            plt.show()

        return ax


