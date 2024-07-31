# Special implementation for 'x_p'

import matplotlib.pyplot as plt
import numpy as np

from . ScalarQuantity import ScalarQuantity
from . OutputException import OutputException

class SPIShardPositions(ScalarQuantity):


    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)
        

    def plotRadialCoordinate(self, shards=None,**kwargs):
        """ 
        Wrapper for ScalarQuantity.plot(), calculating 
        the radial coordinate of the shards instead of 
        the cartesian coordinates. Also allows the user 
        to choose which shards whose radial coordinates 
        should be plotted. 
        NOTE: Currently only valid for cylindrical geometry!
        
        :param slice shards: Shards wose radii should be plotted
        
        :return: Axis object containing the plot
        """
        rhop, _ = self.calcRadialCoordinate(shards)
                
        _rhop = ScalarQuantity(name='\\rho_p', data=rhop, grid=self.grid, output=self.output)

        return _rhop.plot(**kwargs)
        

    def calcRadialCoordinate(self, shards=None, t=None):
        """ 
        Calculates the radial coordinates of the shards 
        (instead of the cartesian coordinates)
        
        :param slice shards: Shards wose radial coordinates should be calculated
        :param slice t: time steps at which the radial coordinates should be calculated
        
        :return: radial coordinate ``rhop`` and poloidal angle ``thetap``.
        """
        if shards is None:
            shards = slice(None)
            
        if t is None:
            t = slice(None)
        
        xp = self.data[:,0::3,0] 
        yp = self.data[:,1::3,0]
        zp = self.data[:,2::3,0]
        
        rhop   = np.sqrt(xp[t,shards]**2+yp[t,shards]**2)
        thetap = np.arctan2(yp[t,shards],xp[t,shards])

        return rhop, thetap


    def plotTrajectoryPoloidal(self, shards=None, ax=None, show=None, color=None):
        """
        Plot the trajectory of one or more shards in a poloidal cross-section.
        """
        black = (87/255, 117/255, 144/255)
        red = (249/255, 65/255, 68/255)

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        if 'eq' not in self.grid:
            raise OutputException("Cannot plot poloidal trajectory when equilibrium data is not stored in output.")

        eq = self.grid.eq

        xp = self.data[:,0::3,0]
        yp = self.data[:,1::3,0]

        if shards is None:
            i1 = np.argmin(yp[-1,:])
            i2 = np.argmax(yp[-1,:])
            shards = (i1, i2)

            if color is None:
                color = red

        if np.isinf(eq.R0): R0 = 1
        else: R0 = eq.R0

        #ax.plot(R0*eq.ROverR0_f[:,-1] - R0, eq.Z_f[:,-1] - eq.Z0, color=black, linewidth=2)
        eq.visualize(ax=ax, shifted=True, maxis=False)

        ax.plot(xp[0,0], yp[0,0], 'o', color=red)
        if color is None:
            ax.plot(xp[:,shards], yp[:,shards])
        else:
            ax.plot(xp[:,shards], yp[:,shards], color=color)

        if np.isinf(eq.R0):
            ax.set_xlabel('Radius $R-R_0$ (m)')
        else:
            ax.set_xlabel('Major radius $R$ (m)')

        ax.set_ylabel('Height $Z$ (m)')
        ax.axis('equal')

        if show:
            plt.show()

        return ax
        
        
    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of shards) covered by this
        quantity. The total number of elements in 'self.data' is the size of
        the grid on which this quantity lives (i.e. 1, since the shards are
        defined on a scalar grid) times this number.
        """
        return self.data.shape[1]
        

