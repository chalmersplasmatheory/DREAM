# Special implementation for 'x_p'

from matplotlib import path
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
        

    def arrivalTime(self, shard=None):
        """
        Estimate the time at which the pellet arrives to the plasma edge.
        """
        if 'eq' not in self.grid:
            raise OutputException("Cannot plot poloidal trajectory when equilibrium data is not stored in output.")

        RMinusR0 = self.grid.eq.RMinusR0_f
        Z = self.grid.eq.ZMinusZ0_f
        ntheta = self.grid.eq.theta.size

        vertices = [(RMinusR0[i,-1], Z[i,-1]) for i in range(ntheta)]
        p = path.Path(vertices)

        xp = self.data[:,0::3,0]
        yp = self.data[:,1::3,0]

        # Check if the shards start within the plasma
        for i in range(xp.shape[1]):
            if p.contains_point((xp[0,i], yp[0,i])):
                return 0

        # Find the pellet which is travelling the fastest in the x direction
        if shard is None:
            shard = np.argmax(np.abs(xp[1,:] - xp[0,:]))
            # Roughly estimate when the fastest pellet reaches the plasma edge
            it_est = np.argmin(np.abs(RMinusR0[0,-1] - xp[:,shard]))
        else:
            xp = xp[:,shard].reshape((xp.shape[0], 1))
            yp = yp[:,shard].reshape((yp.shape[0], 1))
            # Roughly estimate when the fastest pellet reaches the plasma edge
            it_est = np.argmin(np.abs(RMinusR0[0,-1] - xp[:,0]))

        # Determine if the first pellet arrives before or after
        # the estimated time
        arrives_before = False
        for ip in range(xp.shape[1]):
            if p.contains_point((xp[it_est,ip], yp[it_est,ip])):
                arrives_before = True
                break
        
        if arrives_before:
            for it in range(it_est-1, -1, -1):
                inside = False
                for ip in range(xp.shape[1]):
                    if p.contains_point((xp[it,ip], yp[it,ip])):
                        # One pellet has arrived
                        # => break out and go to one time step earlier...
                        inside = True
                        break

                if not inside:
                    return it+1
        else:
            for it in range(it_est+1, self.grid.t.size):
                inside = False
                for ip in range(xp.shape[1]):
                    if p.contains_point((xp[it,ip], yp[it,ip])):
                        inside = True
                        break

                if inside:
                    return it

        return np.nan


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


    def plotAtTime(self, t=-1, shards=None, ax=None, show=None):
        """
        Plot the pellet shards over the poloidal cross-section at a given time.
        """
        black = (87/255, 117/255, 144/255)
        red = (249/255, 65/255, 68/255)

        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        if shards is None:
            shards = slice(None)

        if 'eq' not in self.grid:
            raise OutputException("Cannot plot poloidal trajectory when equilibrium data is not stored in output.")

        eq = self.grid.eq

        xp = self.data[:,0::3,0]
        yp = self.data[:,1::3,0]

        eq.visualize(ax=ax, shifted=True, maxis=False)

        if (xp[0,0] != xp[0,1]) or (yp[0,0] != yp[0,1]):
            print('WARNING: Pellet shards do not start from the same position. Skipping plot of origin.')
        else:
            ax.plot(xp[0,0], yp[0,0], 'o', color=red)

        ax.plot(xp[t,shards], yp[t,shards], 'k.')

        ax.set_xlabel('Radius $R-R_0$ (m)')
        ax.set_ylabel('Height $Z-Z_0$ (m)')
        ax.axis('equal')

        if show:
            plt.show()

        return ax


    def plotTrajectoryPoloidal(self, shards=None, ax=None, show=None, color=None):
        """
        Plot the trajectory of one or more shards in a poloidal cross-section.
        """
        black = (87/255, 117/255, 144/255)
        red = (249/255, 65/255, 68/255)

        if ax is None:
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

        eq.visualize(ax=ax, shifted=True, maxis=False)

        ax.plot(xp[0,0], yp[0,0], 'o', color=red)
        if color is None:
            ax.plot(xp[:,shards], yp[:,shards])
        else:
            ax.plot(xp[:,shards], yp[:,shards], color=color)

        ax.set_xlabel('Radius $R-R_0$ (m)')
        ax.set_ylabel('Height $Z-Z_0$ (m)')
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
        

