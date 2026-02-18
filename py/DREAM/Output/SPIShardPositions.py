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


    def plotAtTime(self, t=-1, shards=None, ax=None, show=None, backgroundQuantity = None, displayGrid = True, displayDrift = False, shardColor = 'k', driftColor = 'k', depColor = 'k', scaleShardsWithSize = False, sizeFactor=5e3, **kwargs):
        """
        Plot the pellet shards over the poloidal cross-section at a given time.
        
        :param int t: Time index to plot.
        :param slice shards: Shards which should be plotted.
        :param matplotlib.pyplot.axis ax: Matplotlib axes object to use for plotting.
        :param bool show: If 'True', shows the plot immediately via a call to
            'matplotlib.pyplot.show()' with 'block=False'. If 'None', this is 
            interpreted as 'True' if 'ax' is also 'None'.
        :param DREAM.Output.FluidQuantity.FluidQuantity backgroundQuantity: 
            FluidQuantity object for which the poloidal contours should be 
            included in the background.
        :param bool displayGrid: Specify wether or not to display the grid cells 
            in the plot.
        :param bool displayGrid: Specify wether or not to display lines illustrating
            the plasmoid drifts in the plot.
        :param shardColor: Color of the plotted shards.
        :param driftColor: Color of the lines illustrating the plasmoid drifts.
        :param depColor: Color of the points marking the deposition position.
        :param bool scaleShardsWithSize: Specify wether or not to scale the marker 
            size of the shards with their size in the simulation.
        :param float sizeFactor: factor used to scale up the shard radii (in meters) 
            to make them visible in the plot
            
        :return: Axis object containing the plot
        """
        black = (87/255, 117/255, 144/255)
        red = (249/255, 65/255, 68/255)

        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        # Plot color scale showing the chosen background quantity, if any
        if backgroundQuantity is not None:
            # We set zorder = 0 to make sure the background color scale is actually plotted in the background and does not cover the shards
            contours, cb = backgroundQuantity.plotPoloidal(ax=ax,show=False, t=t, shifted = True, zorder = 0, **kwargs)
        else:
            contours = None
            
        if shards is None:
            shards = slice(None)

        if 'eq' not in self.grid:
            raise OutputException("Cannot plot poloidal trajectory when equilibrium data is not stored in output.")

        eq = self.grid.eq

        # Retrieve shard position data in cartesian SPI coordinates
        xp = self.data[:,0::3,0]
        yp = self.data[:,1::3,0]
        zp = self.data[:,2::3,0]
        
        # Calculate cylindrical RZ coordinates for the shards
        Rp = np.sqrt((xp+eq.R0)**2+zp**2)
        Zp = yp + eq.Z0

        # If the drift should be displayed, calculate the major radius coordinates of the shifted deposition locations
        if displayDrift:
            Rp_dep = Rp[1:,:] + self.output.other.scalar.ablationDriftMajorRadius.data[:,:,0]

            # The drift is not calculated in the zeroth time steps, so we duplicate the first time step 
            # to make the size of Rp_dep the same as Rp and Zp
            Rp_dep = np.vstack((Rp_dep[0,:].reshape(1,-1),Rp_dep))

        # Display grid if requested
        if displayGrid:
            eq.visualize(ax=ax, shifted=True, maxis=False)
        else:
            # Even if the grid is not displayed, we stil display the outermost flux surface
            ax.plot(eq.RMinusR0_f[:,-1], eq.ZMinusZ0_f[:,-1], color=black, linewidth=2, zorder = 0)
            ax.plot(eq.RMinusR0_f[(0,-1),-1], eq.ZMinusZ0_f[(0,-1),-1], color=black, linewidth=2, zorder = 0) # Close flux surface contour

        # Plot common shard origin, if any
        if (Rp[0,0] != Rp[0,1]) or (Zp[0,0] != Zp[0,1]):
            print('WARNING: Pellet shards do not start from the same position. Skipping plot of origin.')
        else:
            ax.plot(Rp[0,0]-eq.R0, Zp[0,0], 'o', color=red)

        if scaleShardsWithSize:
            # Calculate marker sizes proportional to the shard radii
            rp=self.output.eqsys.Y_p.calcRadii(t=t)
            sizes=rp*sizeFactor

            # Calculate the opacity of the drifts and deposition positions, if requested,
            # such that it increases with the ablation rate
            if displayDrift:

                # Calculate the time derivative Vpdot of the shard volumes, which is proportional to the ablation rate
                if t==0:
                    # If t==0 we plot the first time step, which is the first one where the ablation is calculated
                    Vpdot = 12/5*np.pi*(self.output.eqsys.Y_p.data[t,shards,0]*(self.output.eqsys.Y_p.data[t,shards,0]>=0))*(4/5)*self.output.other.scalar.Ypdot.data[0,shards,0]
                else:
                    # If t>0, we use row t-1 in Ypdot, since Ypdot is shifted compared to Y_p
                    # because Ypdot is not calculated in the zeroth time step
                    Vpdot = 12/5*np.pi*(self.output.eqsys.Y_p.data[t,shards,0]*(self.output.eqsys.Y_p.data[t,shards,0]>=0))*(4/5)*self.output.other.scalar.Ypdot.data[t-1,shards,0]

                # Calculate opacity of the drift and deposition positions
                # This scaling is quite arbitrary chosen, and may not be optimal to quantitatively compare ablation rates,
                # but looks quite good aestetically at least and makes the drift quite clearly visible for most ablating shards
                if np.max(-Vpdot)>0:
                    #alpha = (-Vpdot/np.max(-Vpdot))**0.25
                    alpha = (np.tanh(np.log(-Vpdot/np.max(-Vpdot)*10+1e-30)) + 1)/2
                else:
                    alpha = np.zeros(Vpdot.shape)

        else:
            sizes = 4*np.ones(len(Rp[t,shards]))

            # If the shard marker sizes are not scaled with their actual size, 
            # we also do not vary the opacity of the drifts and deposition locations with the ablation rate
            if displayDrift:
                alpha = np.ones(len(Rp[t,shards]))


        # Plot shards
        ax.scatter(Rp[t,shards]-eq.R0,Zp[t,shards],s=sizes[shards],color=shardColor, edgecolors = 'k', zorder=2)

        # Plot drifts and deposition positions, if requested
        if displayDrift:
            # Cap the drift somewhat outside the plasma
            Rp_dep_max = (eq.RMinusR0_f[0,-1]*1.2+eq.R0)
            Rp_dep[Rp_dep>Rp_dep_max] = Rp_dep_max

            # To make the drifting clouds look more "cloudy", we plot several lines on top of each other
            # with increasing line width, to make the edges less sharp
            for i in range(len(alpha)):
                lw_core = 1.5 # Line width with full opacity
                lw_edge = 3 # Largest line width

                # Number of steps of decreasing opacity (the opacity is the same for all lines individually,
                # but the total opacity increases where several lines are plotted on top of each other)
                nlw = 5

                # Plot lines illustrating the drift.
                # We set zorder = 1 here to plot the drifts behind the shards and deposition locations, to keep them more clearly visible
                for ilw in range(nlw):
                    ax.plot(np.array([Rp[t,shards][i], Rp_dep[t,shards][i]])-eq.R0, np.array([Zp[t,shards][i], Zp[t,shards][i]]), linewidth = lw_core+ilw*(lw_edge-lw_core)/nlw, color = driftColor, alpha = alpha[i]*0.8/nlw, zorder=1, solid_capstyle = 'round')

            # Plot deposition positions. We put these on top of everything else.
            ax.scatter(Rp_dep[t,shards]-eq.R0,Zp[t,shards],s=sizes[shards],color=depColor, alpha = alpha, zorder = 3)


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
        

