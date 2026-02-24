# Base class for kinetic (radius + momentum + time) quantities
#
import warnings

from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np
import scipy

from . OutputException import OutputException
from . UnknownQuantity import UnknownQuantity

from .. import GeriMap
from .. Settings.MomentumGrid import TYPE_PXI
from . PXiGrid import PXiGrid
from . PparPperpGrid import PparPperpGrid

from .. import helpers

class KineticQuantity(UnknownQuantity):

    def __init__(self, name, data, grid, output, momentumgrid=None, attr=list()):
        """
        Constructor.
        """
        super(KineticQuantity, self).__init__(name=name, data=data, attr=attr, grid=grid, output=output)

        self.momentumgrid = momentumgrid

        if momentumgrid is None:
            warnings.warn(
                "KineticQuantity generated without momentumgrid "
                "(e.g. by operating on two kinetic quantities). "
                "This object will have limited functionality."
            )
        else:
            # Cell or flux grid?
            if momentumgrid.p1.size == data.shape[3]:
                self.p1 = momentumgrid.p1
            elif momentumgrid.p1_f.size == data.shape[3]:
                self.p1 = momentumgrid.p1_f
            else:
                raise Exception("Unrecognized shape of data: {}. Expected (nt, nr, np2, np1) = ({}, {}, {}, {}).".format(data.shape, grid.t.size, grid.r.size, momentumgrid.p2.size, momentumgrid.p1.size))

            if momentumgrid.p2.size == data.shape[2]:
                self.p2 = momentumgrid.p2
            elif momentumgrid.p2_f.size == data.shape[2]:
                self.p2 = momentumgrid.p2_f
            else:
                raise Exception("Unrecognized shape of data: {}. Expected (nt, nr, np2, np1) = ({}, {}, {}, {}).".format(data.shape, grid.t.size, grid.r.size, momentumgrid.p2.size, momentumgrid.p1.size))

        if grid.r.size == data.shape[1]:
            self.radius = grid.r
        elif grid.r_f.size == data.shape[1]:
            self.radius = grid.r_f
        else:
            raise Exception("Unrecognized shape of data: {}. Expected (nt, nr, np2, np1) = ({}, {}, {}, {}).".format(data.shape, grid.t.size, grid.r.size, momentumgrid.p2.size, momentumgrid.p1.size))

        self.time = self.grid.t


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        s = f'({self.name}) Kinetic quantity of size NT x NR x NP2 x NP1 = {self.data.shape[0]} x {self.data.shape[1]} x {self.data.shape[2]} x {self.data.shape[3]}'

        if hasattr(self, 'description'):
            s += f'\n:: {self.description}'
        if hasattr(self, 'description_eqn'):
            s += f'\n:: Evolved using: {self.description_eqn}'

        s += f'\n{self.data}'

        return s


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None, momentumgrid=None):
        """
        Creates a new object of the same type where the provided quantities replace
        those of self.
        """
        if name is None:
            name = self.name
        if data is None:
            data = self.data
        if grid is None:
            grid = self.grid
        if output is None:
            output = self.output
        if attr is None:
            attr = {}
            if hasattr(self, "description"):
                attr["description"] = self.description
            if hasattr(self, "description_eqn"):
                attr["equation"] = self.description_eqn
        if momentumgrid is None:
            momentumgrid = self.momentumgrid

        return type(self)(name=name, data=data, grid=grid, output=output, momentumgrid=momentumgrid, attr=attr)

    def angleAveraged(self, t=None, r=None, moment='distribution'):
        r"""
        Returns the angle-averaged distribution function. Depending on
        the input parameters, the whole or only some parts of the spatiotemporal
        distribution can be angle-averaged.

        This method can only be applied to distributions defined on p/xi
        momentum grids.

        Supported moments:
        
        - ``distribution``: :math:`\left\langle f \right\rangle_{\xi_0}`
        - ``density``: :math:`\left\langle V'f\right\rangle_{\xi_0}`
        - ``current``: :math:`\left\langle v\xi_0 V' f\right\rangle_{\xi_0}`
        - ...or a vector (or scalar) to weight the distribution function with.

        where :math:`\left\langle X \right\rangle_{\xi_0} = \int_{-1}^1 X\,\mathrm{d}\xi_0`.
        """
        if self.momentumgrid is None or self.momentumgrid.type != TYPE_PXI:
            raise OutputException("The angle average can only be calculated on p/xi grids.")

        if t is None: t = slice(None)
        if r is None: r = slice(None)
        
        data = self.data[t,r,:]

        if type(moment) == str:
            if moment == 'distribution':
                # Divide by range of xi
                data = data/2
            elif moment == 'density':
                data = data * self.momentumgrid.Vprime_VpVol[r,:]
            elif moment == 'current':
                vPar = self.momentumgrid.getBounceAveragedVpar()
                data = data * vPar[r,:] * self.momentumgrid.Vprime_VpVol[r,:] * scipy.constants.e
        elif type(moment) == float or type(moment) == np.ndarray:
            data = data * moment * self.momentumgrid.Vprime_VpVol[r,:]
        else:
            raise OutputException("Invalid type of parameter 'moment'.")
            
        favg = np.sum(data * self.momentumgrid.DP2[r,:], axis=data.ndim-2)

        return favg

    
    def animate(self, keep=[], r=0, ax=None, repeat=False, repeat_delay=None, speed=None, blit=True, moment='distribution', save=None, dpi=None, **kwargs):
        """
        Creates an animation of the time evolution of the angle average of
        this kinetic quantity.

        :param list keep:        List of time indices to keep in the plot.
        :param r:                Radius to plot angle average for.
        :param ax:               Axes object to use for plotting.
        :param bool repeat:      If ``True``, repeats animation after it is finished.
        :param int repeat_delay: Number of milliseconds to wait before repeating animation.
        :param int speed:        Number of milliseconds to show each frame.
        :param bool blit:        If ``True``, use blitting to optimize drawing.
        :param str moment:       Moment of distribution function to plot (same values as for :py:meth:`DREAM.Output.KineticQuantity.KineticQuantity.angleAveraged`).
        :param str save:         If provided, saves the animation to the named file instead of showing it.
        :param int dpi:          Video resolution (if saving animation to file).
        :param kwargs:           Keyword arguments passed to ``ax.plot()``.
        """
        show = ax is None

        fig = None
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        favg = self.angleAveraged(r=r, moment=moment)
        def update_ani(num, kq, ax, lines, lbl, favg, keeplines, tfac, tunit, keep):
            lbl.set_text(r't = {:.3f} {}'.format(kq.time[num]*tfac, tunit))

            if keep is not None and num in keep:
                idx = keep.index(num)
            else:
                idx = None

            # Iterate over radii
            n = len(lines)
            for i in range(n):
                if favg.ndim == 3: d = favg[num,i,:]
                else: d = favg[num,:]

                line.set_data(kq.p1, d)

                # Keep line after time step has passed?
                if idx is not None:
                    keeplines[idx*n+i].set_data(kq.p1, d)

            return (line, lbl) + tuple(keeplines)

        if speed is None:
            speed = 50

        # Determine number of radii to plot
        if favg.ndim == 3:
            nr = favg.shape[1]
        elif favg.ndim == 2:
            nr = 1
        else:
            raise OutputException("Invalid number of dimensions selected to animate.")

        # Plot at t=0
        colors = GeriMap.get(N=favg.ndim)
        lines = []
        for i in range(nr):
            # Select data to plot
            if favg.ndim == 3: d = favg[0,i,:]
            else: d = favg[0,:]

            if 'color' not in kwargs:
                line, = ax.semilogy(self.p1, d, color=colors(i/(nr+1)), **kwargs)
            else:
                line, = ax.semilogy(self.p1, d, **kwargs)

            lines.append(line)

        # Create placeholders for the 'keep' lines
        keeplines = []
        if keep is not None:
            for i in range(len(keep)):
                for j in range(nr):
                    if 'color' not in kwargs:
                        l, = ax.plot([], [], linewidth=2, color=colors(j/(nr+1)), **kwargs)
                    else:
                        l, = ax.plot([], [], linewidth=2, **kwargs)

                    keeplines.append(l)

        # Set x/y limits
        fmax = np.amax(favg)
        xmin, xmax = self.p1[0], self.p1[-1]
        ymin, ymax = 1e-30*fmax, 10*fmax
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin ,ymax])

        # Determine the relevant time scale
        tmax = self.time[-1]
        idx  = 0
        tfac = 1
        tunits = ['s', 'ms', 'µs', 'ns', 'ps']
        while tmax*tfac < 1 and idx < len(tunits)-1:
            idx += 1
            tfac = (1e3)**(idx)

        xp, yp = 0.60, 0.93
        lymin, lymax = np.log10(ymin), np.log10(ymax)
        tx = xmin+xp*(xmax-xmin)
        ty = lymin+yp*(lymax-lymin)
        txt = ax.text(tx, 10**ty, r't = {:.3f} {}'.format(self.time[0]*tfac, tunits[idx]), usetex=False)

        ax.set_xlabel(r'$r/a$ (m)')

        # Create the animation
        ani = animation.FuncAnimation(fig, update_ani, frames=self.time.size,
            interval=speed, repeat_delay=repeat_delay, repeat=repeat, blit=blit,
            fargs=(self, ax, lines, txt, favg, keeplines, tfac, tunits[idx], keep))

        # Save animation?
        if save:
            writer = animation.FFMpegFileWriter(fps=1000/speed)
            writer.setup(fig, save, dpi=dpi)
            ani.save(save, writer=writer)
            print("Done saving video to '{}'.".format(save))

        if show:
            plt.show()

        return ani


    def get(self, t=None, r=None, p2=None, p1=None):
        """
        Returns data using the specified indexing. If an argument is ``None``,
        this method will return all elements along that dimension.
        """
        sel = [slice(None)] * 4

        if t  is not None: sel[0] = t
        if r  is not None: sel[1] = r
        if p2 is not None: sel[2] = p2
        if p1 is not None: sel[3] = p1

        return self.data[tuple(sel)]


    def moment(self, weight, t=None, r=None):
        """
        Evaluate a moment of this distribution function with the given weighting
        factor.
        """
        if t is None:
            t = range(len(self.time))
        if r is None:
            r = range(len(self.grid.r))

        if np.isscalar(t):
            t = np.asarray([t])
        if np.isscalar(r):
            r = np.asarray([r])

        if np.ndim(weight) != 4:
            _weight = np.ones((self.time.size,self.grid.r.size,self.momentumgrid.p2.size,self.momentumgrid.p1.size))
            
            if np.ndim(weight) == 0:
                weight = _weight*weight
            if np.ndim(weight) == 1:
                weight = _weight*weight[np.newaxis,np.newaxis,np.newaxis,:]
            elif np.ndim(weight) == 2:
                weight = _weight*weight[np.newaxis,np.newaxis,:]
            elif np.ndim(weight) == 3:
                weight = _weight*weight[np.newaxis,:]
                
        q = np.zeros((len(t), len(r)))
        for iT in range(len(t)):
            for iR in range(len(r)):
                q[iT,iR] = self.momentumgrid.integrate2D(self.data[t[iT],r[iR],:] * weight[t[iT],r[iR],:])[iR]
        
        return q

    def plotMesh(self, t=-1, r=0, ax=None, show=None, logarithmic=False, coordinates=None, phaseSpaceWeight=False, **kwargs):
        """
        Visualize this kinetic quantity at one time and radius using a pcolormesh plot,
        representing the cell values as used inside the DREAM solver.

        :param t:           Time index to visualize quantity at.
        :param r:           Radial index to visualize quantity at.
        :param ax:          Matplotlib Axes object to draw plot on.
        :param show:        If ``True``, calls ``matplotlib.pyplot.show()`` with ``block=False`` after plotting the quantity.
        :param logarithmic: If ``True``, plots the base-10 logarithm of the quantity. Default: False.
        :param coordinates: Determines which momentum coordinate system to use. Supported options: "spherical" or "cylindrical". If None (default), chooses the coordinate system used inside DREAM.
        :param phaseSpaceWeight: If ``True``, adds the phase space weight so that the "area in the plot" represents the density. Default: False.
        :param kwargs:      Keyword arguments passed on to ``matplotlib.Axes.pcolormesh()``.
        """
        mg = self.momentumgrid

        if self.momentumgrid is None:
            raise OutputException(
                "Unable to plot kinetic quantity as its momentum grid has not been specified."
            )

        genax = ax is None
        if genax:
            ax = plt.axes()
            if show is None:
                show = True

        data = self.data[t, r, :]
        if data.ndim > 2:
            raise OutputException(
                "Data dimensionality is too high. Unable to visualize kinetic quantity."
            )
        elif data.ndim < 2:
            raise OutputException(
                "Data dimensionality is too low. Unable to visualize kinetic quantity."
            )
        weight_label = ""
        if phaseSpaceWeight:
            weight_label = "(Vp/Vpvol)*"

        if coordinates is None:
            P1, P2 = np.meshgrid(mg.p1_f, mg.p2_f)
            xlabel = mg.p1TeXname
            ylabel = mg.p2TeXname
            if phaseSpaceWeight:
                data = data * mg.Vprime_VpVol[r]
        elif coordinates.lower() == "cylindrical"[: len(coordinates)]:
            P1, P2 = mg.PPAR_f, mg.PPERP_f
            xlabel = PparPperpGrid.P1_TEX_NAME
            ylabel = PparPperpGrid.P2_TEX_NAME
            if phaseSpaceWeight:
                data = data * mg.VprimeCylindrical[r]
        elif coordinates.lower() == "spherical"[: len(coordinates)]:
            P1, P2 = mg.P_f, mg.XI_f
            xlabel = PXiGrid.P1_TEX_NAME
            ylabel = PXiGrid.P2_TEX_NAME
            if phaseSpaceWeight:
                data = data * mg.VprimeSpherical[r]
        else:
            raise ValueError(f"Unrecognized coordinates: {coordinates}")

        sign = ""
        if logarithmic:
            if np.all(data <= 0):
                sign = "$-$"
                data = np.log10(-data)
            else:
                data = np.log10(data)

        pcm = ax.pcolormesh(P1, P2, data, cmap="GeriMap", **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{sign}{weight_label}{self.getTeXName()}")

        if genax:
            plt.colorbar(mappable=pcm, ax=ax)

        if show:
            plt.show(block=False)

        return ax, pcm

    def plot(self, t=-1, r=0, ax=None, show=None, logarithmic=False, coordinates=None, interpolateCylindrical=False, phaseSpaceWeight=False, **kwargs):
        """
        Visualize this kinetic quantity at one time and radius using a filled
        contour plot.

        :param t:           Time index to visualize quantity at.
        :param r:           Radial index to visualize quantity at.
        :param ax:          Matplotlib Axes object to draw plot on.
        :param show:        If ``True``, calls ``matplotlib.pyplot.show()`` with ``block=False`` after plotting the quantity.
        :param logarithmic: If ``True``, plots the base-10 logarithm of the quantity. Default: False.
        :param coordinates: Determines which momentum coordinate system to use. Supported options: "spherical" or "cylindrical". If None (default), chooses the coordinate system used inside DREAM.
        :param interpolateCylindrical: If True and coordinates=="cylindrical", will extend data to pperp=0 to avoid the gaps that form since there's otherwise no cell center there.
        :param phaseSpaceWeight: If ``True``, adds the phase space weight so that the "area in the plot" represents the electron density. Default: False.
        :param kwargs:      Keyword arguments passed on to ``matplotlib.Axes.contourf()``.
        """
        if self.momentumgrid is None:
            raise OutputException("Unable to plot kinetic quantity as its momentum grid has not been specified.")

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        weight_label = ""
        data = self.data[t, r, :]
        if phaseSpaceWeight:
            weight_label = "(Vp/Vpvol)*"
            if coordinates == 'cylindrical'[:len(coordinates)]:
                data = data * self.momentumgrid.VprimeCylindrical[r]
            else:
                data = data * self.momentumgrid.VprimeSpherical[r]

        sign = ''
        if logarithmic:
            if np.all(data <=0):
                sign = '$-$'
                data = np.log10(-data)
            else:
                data = np.log10(data)

        if data.ndim != 2:
            raise OutputException("Data dimensionality is too high. Unable to visualize kinetic quantity.")

        if coordinates is None:
            cp = ax.contourf(self.p1, self.p2, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(self.momentumgrid.getP1TeXName())
            ax.set_ylabel(self.momentumgrid.getP2TeXName())
        # Accept 'spherical' or 'spherica' or 'spheric' or ... 's':
        elif coordinates == 'spherical'[:len(coordinates)]:
            cp = ax.contourf(self.momentumgrid.P, self.momentumgrid.XI, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'$p/mc$')
            ax.set_ylabel(r'$\xi$')
        elif coordinates == 'cylindrical'[:len(coordinates)]:
            if data.shape[1] == self.momentumgrid.PPAR.shape[1] + 1:
                data = (data[:,1:] + data[:,:-1]) / 2
            if data.shape[0] == self.momentumgrid.PPAR.shape[0] + 1:
                data = (data[1:,:] + data[:-1, :]) / 2

            if interpolateCylindrical:
                pperp = np.concatenate((np.zeros(self.momentumgrid.PPERP.shape[1]).reshape((1,-1)), self.momentumgrid.PPERP, np.zeros(self.momentumgrid.PPERP.shape[1]).reshape((1,-1))), axis=0)
                ppar_trailing = (self.momentumgrid.PPAR[-1,:] + (self.momentumgrid.PPAR[-1,:] - self.momentumgrid.PPAR[-2,:])/2).reshape((1,-1))
                ppar_leading = (self.momentumgrid.PPAR[0,:] + (self.momentumgrid.PPAR[0,:] - self.momentumgrid.PPAR[1,:])/2).reshape((1,-1))
                ppar = np.concatenate((ppar_leading, self.momentumgrid.PPAR, ppar_trailing), axis=0)
                
                data_int = np.concatenate((data[0,:].reshape((1,-1)), data, data[-1,:].reshape((1,-1))), axis=0)
                
                cp = ax.contourf(ppar, pperp, data_int, cmap='GeriMap', **kwargs)
            else:
                cp = ax.contourf(self.momentumgrid.PPAR, self.momentumgrid.PPERP, data, cmap='GeriMap', **kwargs)
            ax.set_xlabel(r'$p_\parallel/mc$')
            ax.set_ylabel(r'$p_\perp/mc$')
        else:
            raise OutputException("Unrecognized coordinate type: '{}'.".format(coordinates))

        ax.set_title(f'{sign}{weight_label}{self.getTeXName()}')

        if genax:
            plt.colorbar(mappable=cp, ax=ax)

        if show:
            plt.show(block=False)

        return ax, cp
        
        
        
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
            
        
    def animatePolar(self, t=None, repeat=False, repeat_delay=None, speed=None, dpi=100, save=None,**kwargs):
        """
        Make an animation of poloidal plots of the present quantity, 
        including the specified time steps.
        

        :param slice t: time steps to include in the animation
        :param bool repeat: If ``True``, repeats the animation.
        :param int repeat_delay: Time between consecutive animation runs in milliseconds
        :param int speed: delay between frames in milliseconds
        :param int dpi: animation resolution
        :param str save: title of the file (if any) into which the animation is saved
        """
        
        fig, ax=plt.subplots(1,1)
        
        if t is None:
            t=range(len(self.grid.t))
            
        ax,cb=self.plotPolar(show=False,t=0,**kwargs)
        
        def update_ani(t, kq, ax):
            ax.clear()
            ax=kq.plotPolar(colorbar=False, show=False, t=t,**kwargs)
            
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
		

    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of ion species and
        charge states) covered by this quantity. The total number of elements
        in 'self.data' is the size of the grid on which this quantity lives
        (i.e. scalar grid, fluid grid, or a kinetic grid) times this number.
        """
        return 1

    def resolutionSensors(self, axis, y_floor=0.0, quantityWeighted=False):
        """
        Returns two KineticQuantities representing resolution sensors S_grad and S_curv.

            S_grad: dx * |dy/dx| / (|y| + y_floor)
            S_curv: dx^2 * |d^2y/dx^2| / (|y| + y_floor)

        where y is this quantity's data, and x is the coordinate along `axis`.
        y_floor has the effect of suppressing the sensors where |y| << y_floor.

        S_grad represents the relative amount that data changes between adjacent cells in the `axis` direction.
        S_curv is a similar measure in the second derivative (roughly the relative non-linear change between cells).

        If quantityWeighted is True, returns instead |y| * S_grad, |y| * S_curv,
        representing roughly the absolute change between adjacent cells.
        """
        if self.data.ndim != 4:
            raise OutputException(
                "Resolution sensor assumes 4D data. Found {} for quantity {}".format(
                    self.data.ndim, self.name
                )
            )
        mg = self.momentumgrid
        if axis == 0:
            raise OutputException(
                "Cannot evaluate resolution sensitivity w.r.t. time for quantity {}.".format(
                    self.name
                )
            )
        elif axis == 1:
            x = mg.r
            dx = mg.dr
            xname = "r"
        elif axis == 2:
            xname = mg.p2name
            x = mg.p2
            dx = mg.dp2
        elif axis == 3:
            xname = mg.p1name
            x = mg.p1
            dx = mg.dp1
        else:
            raise ValueError(f"Unrecognized 'axis': {axis}")

        S_grad, S_curv = helpers.compute_S_grad_S_curv(
            self.data, x, dx, axis=axis, y_floor=y_floor
        )

        rel_label = "Relative "
        if quantityWeighted:
            rel_label = ""
            S_grad = S_grad * np.abs(self.data)
            S_curv = S_curv * np.abs(self.data)

        S_grad_qty = KineticQuantity(
            "$S_\\mathrm{{grad}}^{{{}}}$({})".format(xname, self.name),
            data=S_grad,
            grid=self.grid,
            output=self.output,
            momentumgrid=self.momentumgrid,
            attr={
                "description": "{}{}-direction cell gradient uncertainty measure of {}".format(
                    rel_label, xname, self.name
                )
            },
        )
        S_curv_qty = KineticQuantity(
            "$S_\\mathrm{{curv}}^{{{}}}$({})".format(xname, self.name),
            data=S_curv,
            grid=self.grid,
            output=self.output,
            momentumgrid=self.momentumgrid,
            attr={
                "description": "{}{}-direction cell hessian uncertainty measure of {}".format(
                    rel_label, xname, self.name
                )
            },
        )
        return S_grad_qty, S_curv_qty
