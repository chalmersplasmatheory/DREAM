# Base class for fluid (radius + time) quantities
#

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from . OutputException import OutputException
from . UnknownQuantity import UnknownQuantity


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


    def animate(self, keep=[], ax=None, repeat=False, repeat_delay=None, interval=None, blit=True, **kwargs):
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
        if interval is None:
            interval = 50
        
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
        ymin, ymax = 1.1*np.amin(self.data), 1.1*np.amax(self.data)
        if ymin >= 0:
            ymin, ymax = 0, 1.1*np.amax(self.data)

        ax.set_ylim([ymin, ymax])
        
        # Determine relevant time scale
        tmax = self.time[-1]
        idx  = 0
        tfac = 1
        tunits = ['s', 'ms', 'µs', 'ns', 'ps']
        while tmax*tfac < 1 and idx < len(tunits):
            idx += 1
            tfac = (1e3)**(idx)

        xp, yp = 0.03, 0.93
        txt = ax.text(xmin+xp*(xmax-xmin), ymin+yp*(ymax-ymin), r't = {:.3f} {}'.format(self.time[0]*tfac, tunits[idx]), usetex=False)

        ax.set_xlabel(r'$r/a$ (m)')

        # Create the animation
        ani = animation.FuncAnimation(fig, update_ani, frames=self.time.size,
            interval=interval, repeat_delay=repeat_delay, repeat=repeat, blit=blit,
            fargs=(self, ax, line, txt, keeplines, tfac, tunits[idx], keep))

        if show:
            plt.show()


    def get(self, r=None, t=None):
        """
        Returns the data in the specified time or radial
        point. If neither 'r' nor 't' are given, returns
        the full spatiotemporal evolution of the profile.
        """
        if (r is None) and (t is None):
            return self.data
        elif (r is not None) and (t is None):
            return self.data[:,r]
        elif (r is None) and (t is not None):
            return self.data[t,:]
        else:
            return self.data[t,r]

        
    def plot(self, ax=None, show=None, r=None, t=None, colorbar=True, **kwargs):
        """
        Generate a contour plot of the spatiotemporal evolution
        of this quantity.

        ax:   Matplotlib axes object to use for plotting.
        show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.

        RETURNS a matplotlib axis object and a colorbar object
        (which may be 'None' if not used).
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
            cp = ax.contourf(self.radius, self.time, self.data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'Radius $r$ (m)')
            ax.set_ylabel(r'Time $t$')

            cb = None
            if colorbar:
                cb = plt.colorbar(mappable=cp, ax=ax)

            if show:
                plt.show(block=False)

            return ax, cb
        elif (r is not None) and (t is None):
            return self.plotTimeProfile(r=r, ax=ax, show=show)
        elif (r is None) and (t is not None):
            return self.plotRadialProfile(t=t, ax=ax, show=show)
        else:
            raise OutputException("Cannot plot a scalar value. r = {}, t = {}.".format(r, t))


    def plotRadialProfile(self, t=-1, ax=None, show=None):
        """
        Plot the radial profile of this quantity at the specified
        time slice.

        t: Time index to plot.
        ax:   Matplotlib axes object to use for plotting.
        show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.

        RETURNS a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        if np.isscalar(t):
            t = [t]

        lbls = []
        for it in t:
            ax.plot(self.radius, self.data[it,:])

            # Add legend label
            tval, unit = self.grid.getTimeAndUnit(it)
            lbls.append(r'$t = {:.3f}\,\mathrm{{{}}}$'.format(tval, unit))

        ax.set_xlabel(r'Radius $r$ (m)')
        ax.set_ylabel('{}'.format(self.getTeXName()))

        if len(lbls) > 0:
            ax.legend(lbls)

        if show:
            plt.show(block=False)

        return ax


    def plotTimeProfile(self, r=0, ax=None, show=None):
        """
        Plot the temporal profile of this quantity at the specified
        radius.

        r: Radial index to plot evolution for.
        ax:   Matplotlib axes object to use for plotting.
        show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.

        RETURNS a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        if np.isscalar(r):
            r = [r]

        lbls = []
        for ir in r:
            ax.plot(self.time, self.data[:,ir])

            # Add legend label
            lbls.append(r'$r = {:.3f}\,\mathrm{{m}}$'.format(self.radius[ir]))

        ax.set_xlabel(r'Time $t$')
        ax.set_ylabel('{}'.format(self.getTeXName()))

        if len(lbls) > 1:
            ax.legend(lbls)

        if show:
            plt.show(block=False)

        return ax


    def plotIntegral(self, ax=None, show=None):
        """
        Plot the time evolution of the radial integral of this
        quantity.

        ax:   Matplotlib axes object to use for plotting.
        show: If 'True', shows the plot immediately via a call to
              'matplotlib.pyplot.show()' with 'block=False'. If
              'None', this is interpreted as 'True' if 'ax' is
              also 'None'.

        RETURNS a matplotlib axis object.
        """
        if ax is None:
            ax = plt.axes()

            if show is None:
                show = True

        ax.plot(self.time, self.integral())
        ax.set_xlabel(r'Time $t$')
        ax.set_ylabel('{}'.format(self.getTeXIntegralName()))

        if show:
            plt.show(block=False)

        return ax


    def dumps(self, r=None, t=None):
        return self.get(r=r, t=t).__str__()


    def print(self, r=None, t=None):
        """
        Print the data in this quantity.
        """
        print(self.dumps(r,t))


    def integral(self, t=None, w=1.0):
        """
        Evaluate the volume integral of this fluid quantity
        in the given time step using a trapezoidal rule.

        t: Time step to integrate over. If 'None', integrates
           over radius in every time step. May be a slice.
        w: Weighting function.
        """
        if t is None:
            return self.grid.integrate(self.data, w)
        else:
            return self.grid.integrate(self.data[t,:], w)
        

