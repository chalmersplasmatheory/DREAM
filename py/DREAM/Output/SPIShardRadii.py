# Special implementation for 'r_p'

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import animation

from . ScalarQuantity import ScalarQuantity
from . OutputException import OutputException

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
        
        shards: Shards wose radii should be plotted
        """
        _r_p=ScalarQuantity(name=self.name,data=self.calcRadii(shards), grid=self.grid, output=self.output)
        return _r_p.plot(**kwargs)
        
    def calcRadii(self,shards=None, t=None):
        """
        calculates the actual shard radii (instead of r_p**(5/3) 
        as used in the c++ core)
        
        shards: Shards wose radii should be calculated
        t: time steps at whoch the radii should be calculated
        """
        if shards is None:
            shards=slice(None)
            
        if t is None:
            t=slice(None)
            
        data_rp=(self.data[t,shards,0]*(self.data[t,shards,0]>0))**(3.0/5.0)
        return data_rp.reshape(data_rp.shape[0:2])
        
        
    def calcTotalVolume(self,shards=None, t=None):
        """
        calculates the total volume of the specified shards
        
        shards: Shards whose volume should be calculated
        t: time steps at whoch the radii should be calculated
        """
        if shards is None:
            shards=slice(None)
            
        if t is None:
            t=slice(None)
            
        Vp_tot=np.sum(4*np.pi/3*(self.data[t,shards,0]*(self.data[t,shards,0]>0))**(9.0/5.0),axis=-1)
        return Vp_tot.reshape(-1,1)
        
        
    def plotTotalVolume(self, shards=None, **kwargs):
        """
        Wrapper for ScalarQuantity.plot(), calculating the 
        total volume of the specified shards, 
        
        shards: Shards wose volume should be plotted
        """
        _Vp_tot=ScalarQuantity(name='V_{p,tot}',data=self.calcTotalVolume(shards), grid=self.grid, output=self.output)
        return _Vp_tot.plot(**kwargs)
        
    def plotPoloidal(self, ax=None, show=None, t=-1, displayGrid=False, sizeFactor=5e3, backgroundQuantity=None, **kwargs):
        """
        Plot the position and size of the pellet shards, possibly together with 
        a background poloidal contour plot of another fluid quantity at the 
        specified time step. 
        NOTE: Currently assumes a cylindrical flux surface geometry!
        
        ax:   Matplotlib axes object to use for plotting.
        show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.
        t: Time index to plot
        displayGrid: Specify wether or not to display a polar grid in the plot
        sizeFactor: factor used to scale up the shard radii (in meters) to make 
                    them visible in the plot
        backgroundQuantity: FluidQuantity object for which the poloidal
                    contours should be included in the 
                    background

        RETURNS a matplotlib axis object and a colorbar object
        (which may be 'None' if not used).
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
            _,cb = backgroundQuantity.plotPoloidal(ax=ax,show=False, t=t,displayGrid=displayGrid, **kwargs)
                
                
        data_xp=self.output.eqsys.x_p.data[t,0::3,0]
        data_xp.reshape(data_xp.shape[0:2])
        data_yp=self.output.eqsys.x_p.data[t,1::3,0]
        data_yp.reshape(data_xp.shape[0:2])
        
        data_rp=self.calcRadii(t=t)
        
        rho_p=np.sqrt(data_xp**2+data_yp**2)
        theta_p=np.arctan2(data_yp,data_xp)
        sizes=data_rp*sizeFactor
		
        ax.scatter(theta_p,rho_p,s=sizes,color='c')
		
        if show:
            plt.show(block=False)
            
        return ax, cb

        
        
    def animatePoloidal(self,title='DREAM_animation.mp4', t=None, fps=2, dpi=100,backgroundQuantity=None, **kwargs):
        """
        Make an animation of poloidal plots of the SPI shards, 
        including the specified time steps. It is also possible to
        include poloidal contours of a fluid quantity specified by 
        the backgroudQuantity argument
        
        title: title of the resulting mp4 file
        t: time steps to include in the animation
        fps: frame rate of the animation
        dpi: animation resolution
        backgroundQuantity: FluidQuantity object for which the poloidal
                            contours should be included in the 
                            background
        """
        movie=animation.FFMpegWriter(fps=fps)
        
        fig, ax=plt.subplots(1,1)
        
        if t is None:
            t=range(len(self.grid.t))
            
        
        # Make animation
        with movie.saving(fig,title,dpi):
            for it in t:
            	# Clearing the axis apparently clears also the option to make it poloidal,
            	# so therefore we regenerate the axis at every time step
                ax, cb=self.plotPoloidal(show=False,t=it,backgroundQuantity=backgroundQuantity,**kwargs)
                movie.grab_frame()
                ax.clear()
                if backgroundQuantity is not None:
                    cb.remove()
    
        
        
