# Base class for kinetic (radius + momentum + time) quantities
#

import matplotlib.pyplot as plt
import numpy as np
import scipy

from matplotlib import animation

from . OutputException import OutputException
from . UnknownQuantity import UnknownQuantity

from .. Settings.MomentumGrid import TYPE_PXI, TYPE_PPARPPERP


class KineticQuantity(UnknownQuantity):
    

    def __init__(self, name, data, grid, output, momentumgrid=None, attr=list()):
        """
        Constructor.
        """
        super(KineticQuantity, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output)

        self.momentumgrid = momentumgrid


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Kinetic quantity of size NT x NR x NP2 x NP1 = {} x {} x {} x {}\n:: {}\n:: Evolved using: {}\n{}'.format(self.name, self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3], self.description, self.description_eqn, self.data)


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def angleAveraged(self, t=None, r=None, moment='distribution'):
        """
        Returns the angle-averaged distribution function. Depending on
        the input parameters, the whole or only some parts of the spatiotemporal
        distribution can be angle-averaged.

        This method can only be applied to distributions defined on p/xi
        momentum grids.
        """
        if self.momentumgrid is None or self.momentumgrid.type != TYPE_PXI:
            raise OutputException("The angle average can only be calculated on p/xi grids.")

        if t is None: t = slice(None)
        if r is None: r = slice(None)
        
        data = self.data[t,r,:]

        if type(moment) == str:
            if moment == 'distribution': pass
            elif moment == 'density':
                data = data * self.momentumgrid.Vprime_VpVol[r,:]
            elif moment == 'current':
                vPar = self.momentumgrid.getBounceAveragedVpar()
                data = data * vPar[r,:] * self.momentumgrid.Vprime_VpVol[r,:] * scipy.constants.e
        elif type(moment) == float or type(moment) == np.ndarray:
            data = data * moment * self.momentumgrid.Vprime_VpVol[r,:]
        else:
            raise OutputException("Invalid type of parameter 'moment'.")
            
        favg = np.sum(data * self.momentumgrid.DP2[r,:], axis=data.ndim-2) / np.pi

        return favg


    def get(self, t=None, r=None, p2=None, p1=None):
        """
        Returns data using the specified indexing. If an
        argument is 'None', this method will return all
        elements along that dimension.
        """
        sel = [slice(None)] * 4

        if t  is not None: sel[0] = t
        if r  is not None: sel[1] = r
        if p2 is not None: sel[2] = p2
        if p1 is not None: sel[3] = p1

        return self.data[tuple(sel)]


    def moment(self, weight, t=None, r=None):
        """
        Evaluate a moment of this distribution function with the
        given weighting factor.
        """
        if t is None:
            t = range(len(self.grid.t))
        if r is None:
            r = range(len(self.grid.r))

        if np.isscalar(t):
            t = np.asarray([t])
        if np.isscalar(r):
            r = np.asarray([r])

        q = []
        for iT in range(len(t)):
            qr = []
            for iR in range(len(r)):
                qr.append(self.momentumgrid.integrate2D(self.data[t[iT],r[iR],:] * weight)[0])

            q.append(qr)

        q = np.asarray(q)

        return q


    def plot(self, t=-1, r=0, ax=None, show=None, logarithmic=False, coordinates='spherical', **kwargs):
        """
        Plot this kinetic quantity.
        """
        if self.momentumgrid is None:
            raise OutputException("Unable to plot kinetic quantity as its momentum grid has not been specified.")

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        data = None
        if logarithmic:
            data = np.log10(self.data[t,r,:])
        else:
            data = self.data[t,r,:]

        if data.ndim != 2:
            raise OutputException("Data dimensionality is too high. Unable to visualize kinetic quantity.")

        if coordinates is None:
            cp = ax.contourf(self.momentumgrid.p1, self.momentumgrid.p2, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(self.momentumgrid.getP1TeXName())
            ax.set_ylabel(self.momentumgrid.getP2TeXName())
        # Accept 'spherical' or 'spherica' or 'spheric' or ... 's':
        elif coordinates == 'spherical'[:len(coordinates)]:
            cp = ax.contourf(self.momentumgrid.P, self.momentumgrid.XI, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'$p$')
            ax.set_ylabel(r'$\xi$')
        elif coordinates == 'cylindrical'[:len(coordinates)]:
            cp = ax.contourf(self.momentumgrid.PPAR, self.momentumgrid.PPERP, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'$p_\parallel$')
            ax.set_ylabel(r'$p_\perp$')
        else:
            raise OutputException("Unrecognized coordinate type: '{}'.".format(coordinates))

        cb = None
        if genax:
            cb = plt.colorbar(mappable=cp, ax=ax)

        if show:
            plt.show(block=False)

        return ax
        
        
        
    def plotPolar(self, t=-1, r=0, ax=None, show=None, colorbar=True, displayGrid=False, logarithmic=False, thetaMin=0, thetaMax=np.pi, maxMinScale=True, **kwargs):
        """
        Plot this kinetic quantity on a polar axis.
        
        t: Time index to plot
        ax:   Matplotlib axes object to use for plotting.
        show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.
        colorbar: Specify wether or not to include a colorbar
        displayGrid: Specify wether or not to display a polar grid in the plot
        logarithmic: If 'True', plot logarithm of the data
        thetaMin: Minimum pitch angle included in the plot
        thetaMax: Maximum pitch angle included in the plot
        maxMinScale: If 'True', set tha max and min of the color scale to the 
                     maximum and minimum values of the data stored by this object
                     over all time steps at the radial grid point specified by r

        RETURNS a matplotlib axis object and a colorbar object
        (which may be 'None' if not used).
        """
        if self.momentumgrid is None:
            raise OutputException("Unable to plot kinetic quantity as its momentum grid has not been specified.")

        # As we sometimes do not have very many xi-points, the field line parallel direction can 
        # look "empty" if we just plot for the xi-points in the center of the grid cells
        # Therefore, we extend the pitch-angle grid to cover a whole round (with the first and
        # last point at the same angle) and then extend the data accordingly, so that contourf
        # interpolates over the field line parallel direction
        xi=self.momentumgrid.XI
        pitch_angle=np.concatenate((-np.arccos(xi),np.flip(np.arccos(xi))))
        pitch_angle=np.concatenate([pitch_angle,pitch_angle[:1]])
        
        p=self.momentumgrid.P
        p=np.concatenate((p,p))
        p=np.concatenate([p,p[:1]])

		
        genax = ax is None

        if genax:
            ax = plt.subplot(polar=True)
            ax.set_facecolor('k')
            ax.set_ylim([p[0,0],p[0,-1]])
            ax.set_xlim([thetaMin,thetaMax])
            
            if not displayGrid:
                ax.grid(None)
                ax.set_yticklabels([])
                ax.set_xticklabels([])

            if show is None:
                show = True

        data = None
        if logarithmic:
            data = np.log10(self.data[t,r,:])
        else:
            data = self.data[t,r,:]

        if data.ndim != 2:
            raise OutputException("Data dimensionality is too high. Unable to visualize kinetic quantity.")
            
        # Duplicate data in accordance with the extension of the pitch angle grid
        data_plot=np.concatenate((data,np.flip(data,0)))
        data_plot=np.concatenate((data_plot,data_plot[:1]))
		    
        if maxMinScale:
            if logarithmic:
                cp=ax.contourf(pitch_angle,p,data_plot,cmap='GeriMap',levels=np.linspace(np.log10(np.min(np.abs(self.data[:,r,:]))),np.log10(np.max(np.abs(self.data[:,r,:])))),**kwargs)
            else:
                cp=ax.contourf(pitch_angle,p,data_plot,cmap='GeriMap',levels=np.linspace(np.min(self.data[:,r,:]),np.max(self.data[:,r,:])),**kwargs)
        else:
            cp=ax.contourf(pitch_angle,p,data_plot,cmap='GeriMap',**kwargs)
		   
        cb = None
        if colorbar:
            cb = plt.colorbar(mappable=cp, ax=ax)

        if show:
            plt.show(block=False)

        return ax, cb
            
        
    def animatePolar(self,title='DREAM_animation.mp4', t=None, fps=2, dpi=100, **kwargs):
        """
        Make an animation of poloidal plots of the present quantity, 
        including the specified time steps.
        
        title: title of the resulting mp4 file
        t: time steps to include in the animation
        fps: frame rate of the animation
        dpi: animation resolution
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
                ax,cb=self.plotPolar(show=False,t=it,**kwargs)
                movie.grab_frame()
                ax.clear()
                cb.remove()
		

