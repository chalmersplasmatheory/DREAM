
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


    def setEquilibrium(self, eq):
        """
        Set equilibrium data based on output from DREAM.
        """
        self.R0 = eq['R0']
        self.Z0 = eq['Z0']
        self.ROverR0 = eq['ROverR0']
        self.ROverR0_f = eq['ROverR0_f']
        self.Z = eq['Z']
        self.Z_f = eq['Z_f']
        self.theta = eq['theta']


    def visualize(self, *args, ax=None, show=None, **kwargs):
        """
        Visualize this equilibrium.
        """
        red  = (249/255, 65/255, 68/255)
        black = (87/255, 117/255, 144/255)
        gray = (190/255, 190/255, 190/255)

        # Set up axes (if not already done)
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        R0 = self.R0 if not np.isinf(self.R0) else 1

        # Flux surfaces
        ax.plot(R0*self.ROverR0, self.Z, color=gray, linewidth=1, **kwargs)
        #ax.plot(R0*self.ROverR0_f, self.Z_f, color=gray, linewidth=1, **kwargs)
        # Limiter
        ax.plot(R0*self.ROverR0_f[:,-1], self.Z_f[:,-1], color=black, linewidth=2, **kwargs)
        # Magnetic axis 
        ax.plot(R0, self.Z0, 'x', color=red)
        ax.axis('equal')

        if np.isinf(self.R0):
            ax.set_xlabel('Major radius $R-R_0$ (m)')
        else:
            ax.set_xlabel('Major radius $R$ (m)')
        ax.set_ylabel('Height $Z$ (m)')

        if show:
            plt.show()

