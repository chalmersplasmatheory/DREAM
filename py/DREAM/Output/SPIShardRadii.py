# Special implementation for 'r_p'

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import animation

from . ScalarQuantity import ScalarQuantity
from .FluidQuantity import FluidQuantity
from . OutputException import OutputException


anim_contours = None


class SPIShardRadii(ScalarQuantity):


    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, attr=attr, grid=grid, output=output)
        self.nshard=data.shape[1]
        

    def plotRadii(self, shards=None,**kwargs):
        """
        Wrapper for ScalarQuantity.plot(), calculating the actual 
        shard radii (instead of r_p**(5/3) as used in the c++ core), 
        and also allowing the user to select which shards' radii to plot
        
        :param slice shards: Shards wose radii should be plotted
        
        :return: Axis object containing the plot
        """
        _r_p=ScalarQuantity(name=self.name, data=self.calcRadii(shards), grid=self.grid, output=self.output)
        return _r_p.plot(**kwargs)
        

    def calcRadii(self,shards=None, t=None):
        """
        calculates the actual shard radii (instead of r_p**(5/3) 
        as used in the c++ core)
        
        :param slice shards: Shards wose radii should be calculated
        :param slice t: time steps at whoch the radii should be calculated
        
        :return: shard radii
        """
        if shards is None:
            shards = slice(None)
            
        if t is None:
            t = slice(None)
            
        rp = (self.data[t,shards,0]*(self.data[t,shards,0]>0))**(3.0/5.0)
        return rp
        
        
    def calcTotalVolume(self,shards=None, t=None):
        """
        calculates the total volume of the specified shards
        
        shards: Shards whose volume should be calculated
        t: time steps at which the radii should be calculated
        
        :return: Total volume of all specified shards combined
        """
        if shards is None:
            shards=slice(None)
            
        if t is None:
            t=slice(None)
            
        Vp_tot=np.sum(4*np.pi/3*(self.data[t,shards,0]*(self.data[t,shards,0]>0))**(9.0/5.0),axis=-1)
        return Vp_tot


    def calcTotalParticles(self, shards=None, t=None):
        """
        Calculates the total number of particles in the specified shards.
        """
        return self.calcTotalVolume(shards=shards, t=t) / molarVolume * N_A
        
        
    def plotTotalVolume(self, shards=None, **kwargs):
        """
        Wrapper for ScalarQuantity.plot(), calculating the 
        total volume of the specified shards, 
        
        :param slice shards: Shards wose volume should be plotted
        """
        _Vp_tot=ScalarQuantity(name='V_{p,tot} [m$^3$]',data=self.calcTotalVolume(shards), grid=self.grid, output=self.output)
        return _Vp_tot.plot(**kwargs)

        
    def plotAblatedVolume(self, shards = None, **kwargs):
        nt = len(self.grid.t)
        nr = len(self.grid.r)
        data = np.zeros((nt,nr))
        for it in range(1,nt):
        
            data[it,:] = data[it-1,:]
            
            rpPrev = self.calcRadii(shards = shards, t = it-1)
            rp = self.calcRadii(shards = shards,t = it)
            VPrev = 4*np.pi*rpPrev**3/3
            V = 4*np.pi*rp**3/3
            DV = VPrev - V
            r, theta = self.output.eqsys.x_p.calcRadialCoordinate(shards = shards, t=it) # NOTE: assumes r stays the same during the whole saved time step! Should be improved!
            for ir in range(nr):
                data[it,ir] += np.sum(DV[(r>self.grid.r_f[ir]) & (r<self.grid.r_f[ir+1])])
                
        _ablatedVolume = FluidQuantity(name = "ablated volume [m$^3$]", data=data, attr=list(), grid=self.grid, output=self.output)
        return _ablatedVolume.plot(**kwargs)
            
        
    def plotPoloidal(self, ax=None, show=None, t=-1, displayGrid=False, sizeFactor=5e3, backgroundQuantity=None, return_artists=False, **kwargs):
        """
        Plot the position and size of the pellet shards, possibly together with 
        a background poloidal contour plot of another fluid quantity at the 
        specified time step. 
        NOTE: Currently assumes a cylindrical flux surface geometry!
        
        :param matplotlib.pyplot.axis ax:   Matplotlib axes object to use for plotting.
        :param bool show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.
        :param int t: Time index to plot
        :param bool displayGrid: Specify wether or not to display a polar grid in the plot
        :param float sizeFactor: factor used to scale up the shard radii (in meters) to make 
                    them visible in the plot
        :param DREAM.Output.FluidQuantity.FluidQuantity backgroundQuantity: FluidQuantity object for which the poloidal
                    contours should be included in the 
                    background

        :return: a matplotlib axis object and a colorbar object (which may be 'None' if not used).
        """
        genax = ax is None

        if genax:
            ax = plt.subplot(polar=True)
            ax.set_facecolor('k')
            ax.set_ylim([self.grid.r[0],self.grid.r[-1]])
            ax.set_title('t = '+str(self.grid.t[t]))
            cb=None

            if not displayGrid:
                ax.grid(None)
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                

            if show is None:
                show = True
                
        if backgroundQuantity is not None:
            contours, cb = backgroundQuantity.plotPoloidal(ax=ax,show=False, t=t,displayGrid=displayGrid, **kwargs)
        else:
            contours = None
                
        data_xp=self.output.eqsys.x_p.data[t,0::3,0]
        data_xp.reshape(data_xp.shape[0:2])
        data_yp=self.output.eqsys.x_p.data[t,1::3,0]
        data_yp.reshape(data_xp.shape[0:2])
        
        data_rp=self.calcRadii(t=t)
        
        rho_p, theta_p=self.output.eqsys.x_p.calcRadialCoordinate(t=t)
        sizes=data_rp*sizeFactor
		
        artist = ax.scatter(theta_p,rho_p,s=sizes,color='c')
		
        if show:
            plt.show(block=False)
            
        if return_artists:
            return ax, cb, artist, contours
        else:
            return ax, cb
        
        
    def animatePoloidal(
        self, t=None, repeat=False, repeat_delay=None, speed=None,
        dpi=100, save=None, displayGrid=False, backgroundQuantity=None,
        sizeFactor=5e3, **kwargs):
        """
        Make an animation of poloidal plots of the SPI shards, 
        including the specified time steps. It is also possible to
        include poloidal contours of a fluid quantity specified by 
        the backgroudQuantity argument
        
        :param slice t: time steps to include in the animation
        :param bool repeat: If ``True``, repeats the animation.
        :param int repeat_delay: Time between consecutive animation runs in milliseconds
        :param int speed: delay between frames in milliseconds
        :param bool blit: Whether to use blitting when drawing the frames or not.
        :param float dpi: animation resolution
        :param str save: title of the file (if any) into which the animation is saved
        :param DREAM.Output.FluidQuantity.FluidQuantity backgroundQuantity: FluidQuantity object for which the poloidal contours should be included in the background
        """
        global anim_contours

        fig, ax=plt.subplots(1,1)
        
        if t is None:
            t=range(len(self.grid.t))
            
        ax, cb, shards, anim_contours = self.plotPoloidal(show=False,t=0, backgroundQuantity=backgroundQuantity, return_artists=True, **kwargs)
        
        def update_ani(t, ssr, ax):
            global anim_contours

            if anim_contours is not None:
                for c in anim_contours.collections:
                    c.remove()

                anim_contours, cb = backgroundQuantity.plotPoloidal(ax=ax, show=False, t=t, colorbar=False, displayGrid=displayGrid, **kwargs)

            data_rp = self.calcRadii(t=t)
            sizes = data_rp*sizeFactor
            rho_p, theta_p = self.output.eqsys.x_p.calcRadialCoordinate(t=t)
            shards.set_sizes(sizes)
            shards.set_offsets((theta_p, rho_p))

            return shards, anim_contours

        
        if speed is None:
            speed = 50
            
        # Create the animation
        ani = animation.FuncAnimation(fig, update_ani, frames=t,
            repeat=repeat, repeat_delay=repeat_delay,
            interval=speed, fargs=(self, ax))
        
        if save:
            # Make animation
            writer = animation.FFMpegFileWriter(fps=1000/speed)
            writer.setup(fig, save, dpi=dpi)
            ani.save(save, writer=writer)
            print("Done saving video to '{}'.".format(save))
		            
        plt.show()
    
        
    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of shards) covered by this
        quantity. The total number of elements in 'self.data' is the size of
        the grid on which this quantity lives (i.e. 1, since the shards are
        defined on a scalar grid) times this number.
        """
        return self.data.shape[1]
        

