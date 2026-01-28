# Base class for fluid (radius + time) quantities
#

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import animation

from . OutputException import OutputException
from . UnknownQuantity import UnknownQuantity


anim_contours = None


class FluidQuantity(UnknownQuantity):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super(FluidQuantity, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output)

        # Cell or flux grid?
        if data.shape[1] == self.grid.r.size:
            self.radius = self.grid.r
        elif data.shape[1] == self.grid.r.size+1:
            self.radius = self.grid.r_f
        else:
            raise Exception("Unrecognized shape of data for '{}': {}. Expected (nt, nr) = ({}, {}).".format(name, data.shape, grid.t.size, grid.r.size))

        self.time = self.grid.t

    
    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        s = self.__str__() + "\n"
        if hasattr(self, 'description') and hasattr(self, 'description_eqn'):
            s += ":: {}\n:: Evolved using: {}\n".format(self.description, self.description_eqn)
        s += self.dumps()
        return s


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Fluid quantity of size NT x NR = {} x {}'.format(self.name, self.data.shape[0], self.data.shape[1])


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def animate(self, keep=[], ax=None, repeat=False, repeat_delay=None, speed=None, blit=True, save=None, dpi=None, **kwargs):
        """
        Creates an animation of the time evolution of this
        fluid quantity.

        :param list keep:   List of time indices to keep after plotting.
        :param bool repeat: If ``True``, repeats the animation.
        """
        show = ax is None

        fig = None
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        def update_ani(num, fq, ax, line, lbl, keeplines, tfac, tunit, keep):
            lbl.set_text(r't = {:.3f} {}'.format(fq.time[num]*tfac, tunit))
            line.set_data(fq.radius, self.data[num,:])

            if keep is not None and num in keep:
                idx = keep.index(num)
                keeplines[idx].set_data(fq.radius, self.data[num,:])

            return (line, lbl) + tuple(keeplines)

        # Automatically determine the plotting interval
        if speed is None:
            speed = 50
        
        line, = ax.plot(self.radius, self.data[0,:], 'k', linewidth=2, **kwargs)

        # Create placeholders for the 'keep' lines
        keeplines = []
        if keep is not None:
            for i in range(len(keep)):
                l, = ax.plot([], [], linewidth=2, **kwargs)
                keeplines.append(l)

        xmin, xmax = 0, self.radius[-1]
        ax.set_xlim([xmin, xmax])

        # Set y limits
        data = self.data[:]     # Make sure to load data only once
        ymin, ymax = 1.1*np.amin(data), 1.1*np.amax(data)
        if ymin >= 0:
            ymin, ymax = 0, 1.1*np.amax(data)

        ax.set_ylim([ymin, ymax])
        
        # Determine relevant time scale
        tmax = self.time[-1]
        idx  = 0
        tfac = 1
        tunits = ['s', 'ms', 'µs', 'ns', 'ps']
        while tmax*tfac < 1 and idx < len(tunits)-1:
            idx += 1
            tfac = (1e3)**(idx)

        xp, yp = 0.03, 0.93
        txt = ax.text(xmin+xp*(xmax-xmin), ymin+yp*(ymax-ymin), r't = {:.3f} {}'.format(self.time[0]*tfac, tunits[idx]), usetex=False)

        ax.set_xlabel(r'$r/a$ (m)')

        # Create the animation
        ani = animation.FuncAnimation(fig, update_ani, frames=self.time.size,
            interval=speed, repeat_delay=repeat_delay, repeat=repeat, blit=blit,
            fargs=(self, ax, line, txt, keeplines, tfac, tunits[idx], keep))

        # Save animation?
        if save:
            writer = animation.FFMpegFileWriter(fps=1000/speed)
            writer.setup(fig, save, dpi=dpi)
            ani.save(save, writer=writer)
            print("Done saving video to '{}'.".format(save))

        if show:
            plt.show()


    def get(self, r=None, t=None):
        """
        Returns the data in the specified time or radial point. If neither ``r``
        nor ``t`` are given, returns the full spatiotemporal evolution of the
        profile.
        """
        if (r is None) and (t is None):
            return self.data[:]
        elif (r is not None) and (t is None):
            return self.data[:,r]
        elif (r is None) and (t is not None):
            return self.data[t,:]
        else:
            return self.data[t,r]


    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of ion species and
        charge states) covered by this quantity. The total number of elements
        in 'self.data' is the size of the grid on which this quantity lives
        (i.e. scalar grid, fluid grid, or a kinetic grid) times this number.
        """
        return 1

        
    def plot(self, ax=None, show=None, r=None, t=None, log=False, colorbar=True, VpVol=False, weight=None, unit='s', **kwargs):
        """
        Generate a contour plot of the spatiotemporal evolution of this
        quantity.

        :param ax:       Matplotlib axes object to use for plotting.
        :param show:     If 'True', shows the plot immediately via a call to ``matplotlib.pyplot.show()`` with ``block=False``. If ``None``, this is interpreted as ``True`` if ``ax`` is also ``None``.
        :param log:      If ``True``, plot on a logarithmic scale.
        :param colorbar: If ``True``, and a 2D plot is requested, also draw a colorbar.
        :param VpVol:    Weight quantity with ``grid.VpVol`` when plotting.
        :param weight:   Optional quantity to weight this quantity with when plotting.

        :return: a matplotlib axis object and a colorbar object (which may be 'None' if not used).
        """
        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        # If the data is 1D, make sure it is plotted
        # as such (and not as a contour plot)
        if self.data.shape[0] == 1:
            t = 0
        elif self.data.shape[1] == 1:
            r = 0
        
        if (r is None) and (t is None):
            data = np.copy(self.data[:])
            if VpVol:
                data *= self.grid.VpVol[:]
            if weight is not None:
                data *= weight

            if log:
                data = np.log10(np.abs(data))

            time = self.time * self._getTimeUnitFactor(unit)

            cp = ax.contourf(self.radius, time, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'Radius $r$ (m)')
            ax.set_ylabel(fr'Time $t$ ({unit})')

            cb = None
            if colorbar:
                cb = plt.colorbar(mappable=cp, ax=ax)

            if show:
                plt.show(block=False)

            return ax, cb
        elif (r is not None) and (t is None):
            return self.plotTimeProfile(r=r, ax=ax, show=show, VpVol=VpVol, weight=weight, log=log, **kwargs)
        elif (r is None) and (t is not None):
            return self.plotRadialProfile(t=t, ax=ax, show=show, VpVol=VpVol, weight=weight, log=log, **kwargs)
        else:
            raise OutputException("Cannot plot a scalar value. r = {}, t = {}.".format(r, t))


    def plotPoloidal(self, ax=None, show=None, t=-1, colorbar=True, return_contours=False, displayGrid=False, maxis = False, shifted = True, maxMinScale=True, levels = None, logscale=False, **kwargs):
        """
        Plot the radial profile of this quantity revolved over a
        poloidal cross section at the specified time step.

        :param matplotlib.pyplot.Axis ax:   Matplotlib axes object to use for plotting.
        :param bool show: If 'True', shows the plot immediately via a call to 'matplotlib.pyplot.show()' with 'block=False'. If 'None', this is interpreted as 'True' if 'ax' is also 'None'.
        :param int t: Time index to plot.
        :param matplotlib.pyplot.Colorbar colorbar: Specify wether or not to include a colorbar.
        :param bool return_contours: Specify wether or not to return the contours
        :param bool displayGrid: Specify wether or not to display the flux surfaces separating the grid cells in the plot.
        :param bool maxis: Specify wether or not to display the magnetic axis
        :param bool shifted: If 'True', the origin is shifted to the magnetic axis, otherwise the origin is in the center of the torus
        :param bool maxMinScale: If 'True', set tha max and min of the color scale to the maximum and minimum values of the data stored by this object over all time steps.
        :param numpy.ndarray levels: Levels for the color scale (only used if maxMinScale is False)
        :param bool logscale: If 'True', plot the contours with a logarithmic color scale

        :return: a matplotlib axis object and a colorbar object (which may be 'None' if not used).
        """
        
        black = (87/255, 117/255, 144/255)
        
        genax = ax is None

        if genax:
            ax = plt.axes()
            ax.set_title('t = '+str(self.grid.t[t]))
            if shifted:
                ax.set_xlabel('Radius $R-R_0$ (m)')
                ax.set_ylabel('Height $Z-Z_0$ (m)')
            else:
                ax.set_xlabel('Radius $R$ (m)')
                ax.set_ylabel('Height $Z$ (m)')

            ax.axis('equal')

            if show is None:
                show = True

        if 'eq' not in self.grid:
            raise OutputException("Cannot plot poloidal contrours when equilibrium data is not stored in output.")

        theta = self.grid.eq.theta

        # Shift origin as requested
        if shifted:
            R = self.grid.eq.RMinusR0
            R_f = self.grid.eq.RMinusR0_f
            Z = self.grid.eq.ZMinusZ0
            Z_f = self.grid.eq.ZMinusZ0_f
        else:
            R = self.grid.eq.RMinusR0 + self.grid.eq.R0
            R_f = self.grid.eq.RMinusR0_f + self.grid.eq.R0
            Z = self.grid.eq.ZMinusZ0 + self.grid.eq.Z0
            Z_f = self.grid.eq.ZMinusZ0_f + self.grid.eq.Z0
            
        # Add the innermost and outermost flux surfaces of the flux grid,
        # so that the contours fill the entire plasma cross section
        R = np.hstack((R_f[:,0].reshape(-1,1), R))
        R = np.hstack((R, R_f[:,-1].reshape(-1,1)))
        Z = np.hstack((Z_f[:,0].reshape(-1,1), Z))
        Z = np.hstack((Z, Z_f[:,-1].reshape(-1,1)))
        
        # Create matrix with the data in each radial grid cell copied to all values of theta
        # Take logarithm of data if requested
        if logscale:
        	data_mat=np.log10(self.data[t,:])*np.ones((len(theta),len(self.grid.r)))
        else:
        	data_mat=self.data[t,:]*np.ones((len(theta),len(self.grid.r)))
        	
        # Duplicate the data in the innermost and outermost grid cell centers
        # to be plotted at the edges of the flux grid so that the contours fill the entire plasma cross section
        data_mat = np.hstack((data_mat[:,0].reshape(-1,1),data_mat))
        data_mat = np.hstack((data_mat,data_mat[:,-1].reshape(-1,1)))
        	
    # Calculate contour levels ranging between the maximum and minimum value of this fluidQuantity, if requested
        if maxMinScale:
            if logscale:
                levels=np.linspace(np.min(np.log10(self.data)),np.max(np.log10(self.data)))
            else:
                levels=np.linspace(np.min(self.data),np.max(self.data))

        # Plot contours of this fluidQuantity
        cp = ax.contourf(R, Z, data_mat, cmap='GeriMap', levels = levels, **kwargs)
			
        # Create colorbar, if requested
        cb = None
        if colorbar:
            cb = plt.colorbar(mappable=cp, ax=ax)
            if logscale:
                cb.ax.set_ylabel(r'$\log _{10}($'+'{}'.format(self.getTeXName()+')'))
            else:
                cb.ax.set_ylabel('{}'.format(self.getTeXName()))
                
        # Display grid cells, if requested
        if displayGrid:
            self.grid.eq.visualize(ax=ax, maxis=maxis, shifted = shifted)
        else:
            # We display the edge of the plasma even when we do not display the whole grid
            ax.plot(R[:,-1], Z[:,-1], color=black, linewidth=2, **kwargs)
            ax.plot(R[(0,-1),-1], Z[(0,-1),-1], color=black, linewidth=2, **kwargs) # Close contour
            
        if show:
            plt.show(block=False)
            
        if return_contours:
            return ax, cb, cp
        else:
            return ax, cb

        
    def animatePoloidal(self, t=None, repeat=False, repeat_delay=None, speed=None, dpi=100, save=None,**kwargs):
        """
        Make an animation of poloidal plots of the present quantity, 
        including the specified time steps.
        
        :param slice t: time steps to include in the animation
        :param bool repeat: If ``True``, repeats the animation.
        :param int repeat_delay: Time between consecutive animation runs in milliseconds
        :param int speed: delay between frames in milliseconds
        :param float dpi: animation resolution
        :param str save: title of the file (if any) into which the animation is saved
        """
        global anim_contours
        
        fig, ax=plt.subplots(1,1)
        
        if t is None:
            t=range(len(self.grid.t))
            
        ax, cb, anim_contours = self.plotPoloidal(show=False, t=0, return_contours=True, **kwargs)
        
        def update_ani(t, fq, ax):
            global anim_contours

            for c in anim_contours.collections:
                c.remove()

            anim_contours, cb = fq.plotPoloidal(ax=ax, show=False, t=t, colorbar=False, **kwargs)

            return anim_contours
            
        if speed is None:
            speed = 50

        # Create the animation
        ani = animation.FuncAnimation(fig, update_ani, frames=t,
            repeat=repeat, repeat_delay=repeat_delay, interval=speed,
            fargs=(self, ax))
        
        if save:
            # Make animation
            writer = animation.FFMpegFileWriter(fps=1000/speed)
            writer.setup(fig, save, dpi=dpi)
            ani.save(save, writer=writer)
            print("Done saving video to '{}'.".format(save))
		            
        plt.show()

    def plotRadialProfile(self, t=-1, ax=None, show=None, VpVol=False, weight=None, log=False, **kwargs):
        """
        Plot the radial profile of this quantity at the specified time slice.

        :param t:      Time index to plot.
        :param ax:     Matplotlib axes object to use for plotting.
        :param show:   If ``True``, shows the plot immediately via a call to ``matplotlib.pyplot.show()`` with ``block=False``. If ``None``, this is interpreted as ``True`` if ``ax`` is also ``None``.
        :param VpVol:  If ``True``, weight the radial profile with the spatial jacobian V'.
        :param weight: Optional quantity to weight this quantity with when plotting.
        :param log:    If ``True``, plot on a logarithmic scale.

        :return: a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        if np.isscalar(t):
            t = [t]

        lbls = []
        vpv = self.grid.VpVol[:]
        for it in t:
            data = np.copy(self.data[it,:])
            wlbl = ''
            if VpVol:
                data *= vpv
                wlbl += "*V'"
            if weight is not None:
                data *= weight
                wlbl += '*w'


            if log:
                if np.any(data>0):
                    ax.semilogy(self.radius, data, **kwargs)
                else:
                    ax.semilogy(self.radius, -data, '--', **kwargs)
            else:
                ax.plot(self.radius, data, **kwargs)

            # Add legend label
            tval, unit = self.grid.getTimeAndUnit(it)
            lbls.append(r'$t = {:.3f}\,\mathrm{{{}}}$'.format(tval, unit))

        ax.set_xlabel(r'Radius $r$ (m)')
        ax.set_ylabel('{}{}'.format(self.getTeXName(), wlbl))

        if len(lbls) > 0:
            ax.legend(lbls)

        if show:
            plt.show(block=False)

        return ax   	


    def plotTimeProfile(self, r=0, ax=None, show=None, VpVol=False, weight=None, log=False, **kwargs):
        """
        Plot the temporal profile of this quantity at the specified radius.

        :param r:      Radial index to plot evolution for.
        :param ax:     Matplotlib axes object to use for plotting.
        :param show:   If ``True``, shows the plot immediately via a call to ``matplotlib.pyplot.show()`` with ``block=False``. If ``None``, this is interpreted as ``True`` if ``ax`` is also ``None``.
        :param VpVol:  If ``True``, weight the radial profile with the spatial jacobian V'.
        :param weight: Optional quantity to weight this quantity with when plotting.
        :param log:    If ``True``, plot on a logarithmic scale.

        :return: a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        if np.isscalar(r):
            r = [r]

        lbls = []
        for ir in r:
            data = np.copy(self.data[:,ir])
            wlbl = ''
            if VpVol:
                data *= self.grid.VpVol[ir]
                wlbl += "*V'"
            if weight is not None:
                data *= weight
                wlbl += '*w'

            if log:
                if np.any(data>0):
                    ax.semilogy(self.time, data, **kwargs)
                else:
                    ax.semilogy(self.time, -data, '--', **kwargs)
            else:
                ax.plot(self.time, data, **kwargs)

            # Add legend label
            lbls.append(r'$r = {:.3f}\,\mathrm{{m}}$'.format(self.radius[ir]))

        ax.set_xlabel(r'Time $t$')
        ax.set_ylabel('{}{}'.format(self.getTeXName(), wlbl))

        if len(lbls) > 1:
            ax.legend(lbls)

        if show:
            plt.show(block=False)

        return ax


    def plotIntegral(self, ax=None, show=None, unit='s', time_shift = 0, time_scale_factor = 1.0, w=1.0, time_derivative = False, log=False, **kwargs):
        """
        Plot the time evolution of the radial integral of this quantity.

        :param ax:   Matplotlib axes object to use for plotting.
        :param show: If ``True``, shows the plot immediately via a call to ``matplotlib.pyplot.show()`` with ``block=False``. If ``None``, this is interpreted as ``True`` if ``ax`` is also ``None``.

        :return: a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        time = self.time * self._getTimeUnitFactor(unit)
        time = time + time_shift
        time = time*time_scale_factor

        if time_derivative:
            integrated_data = self.integral(w=w)
            tm = time[1:-1]
            v = (integrated_data[2:]-integrated_data[:-2])/((time[2:]-time[:-2])/time_scale_factor)
        else:
            tm = time
            v = self.integral(w=w)

        if log:
            ax.semilogy(tm, v, **kwargs)
        else:
            ax.plot(tm, v, **kwargs)

        ax.set_xlabel(fr'Time $t$ ({unit})')
        ax.set_ylabel('{}'.format(self.getTeXIntegralName()))

        if show:
            plt.show(block=False)

        return ax


    def lineIntegrated(self, t=None, x0=np.array([0,10,0]), n = np.array([0,-1,0]), lmax = 20, nl = 1000, normaliseToPathLength=False):
        """
        Evaluate the line integral of this fluid quantity along a specified line.
        
        :param t: Time steps to evaluate the line integrated density for. If ``None``, the line integral is calculated for all time steps. May be a slice
        :param x0: Starting point for the line to integrate along, in cartesian coordinates centered at the magnetic axis at pho=0 (same coordinates as used to track the SPI shards). One probably wants to set this point outside the plasma (the integrand is set to zero outside the plasma).
        :param n: Vector specifying the direction to integrate along
        :param lmax: Length of the line to integrate over.
        :param nl: number of points along the line used for the numerical integration
        :param normaliseToPathLength: If ``True``, divide the integral by the length of the part of the line residing inside the plasma
        """
        if t is None:
            t=slice(t)

        # define line in Cartesian coordinates (SPI coordinates)
        l = np.linspace(0,lmax,nl)
        dl = l[1]-l[0]
        x = x0[0] + n[0]*l
        y = x0[1] + n[1]*l
        z = x0[2] + n[2]*l

        # get r coordinates for all points along the line
        r, _, _ = self.output.grid.eq.getRThetaPhiFromCartesian(x,y,z)

        # compute line-integral
        lineIntegral = np.zeros((len(self.grid.t[t]),1))
        for i in range(nl):
            if r[i]<self.grid.a:
                lineIntegral+=self.data.data[t, (r[i]>self.grid.r_f[:-1]) & (r[i]<self.grid.r_f[1:])]*dl
        
        # divide by path length?
        if normaliseToPathLength:
            L = np.max(l[r<self.grid.a]) - np.min(l[r<self.grid.a])
            return lineIntegral / L
        
        return lineIntegral


    def dumps(self, r=None, t=None):
        return self.get(r=r, t=t).__str__()


    def print(self, r=None, t=None):
        """
        Print the data in this quantity.
        """
        print(self.dumps(r,t))


    def integral(self, t=None, w=1.0):
        """
        Evaluate the volume integral of this fluid quantity in the given time
        step using a trapezoidal rule.

        :param t: Time step to integrate over. If ``None``, integrates over radius in every time step. May be a slice.
        :param w: Weighting function.
        """
        if t is None:
            return self.grid.integrate(self.data[:], w)
        else:
            return self.grid.integrate(self.data[t,:], w)


    def _getTimeUnitFactor(self, unit):
        """
        Converts a time unit given as a string to a numeric factor
        for converting the 'grid.time' vector to the specified units
        (i.e. from seconds).
        """
        if unit == 's': return 1
        elif unit == 'ms': return 1e3
        elif unit == 'µs': return 1e6
        elif unit == 'ns': return 1e9
        else:
            raise ValueError(f"Unrecognized time unit: '{unit}'.")
        

