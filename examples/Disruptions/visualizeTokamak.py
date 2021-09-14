# Routines for visualizing a tokamak.

from DREAM.DREAMSettings import DREAMSettings
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np


def plotSignal(ax, x, y, name, xlabel=False):
    """
    Plot the given signal nicely.

    :param ax: Axes to plot the signal on.
    :param x:  X coordinates.
    :param y:  Y coordinates.
    """
    ax.plot(x, y, linewidth=2, color=(87/255, 117/255, 144/255))
    if xlabel:
        ax.set_xlabel(r'Minor radius (m)')
    ax.set_ylabel(name)

    ax.set_xlim([x[0], x[-1]])
    ax.set_ylim([0, np.amax(y)*1.2])


def visualize(tokamak, nticks=None, Tticks=None):
    """
    Visualize the given tokamak module.

    :param tokamak: Module reference for the tokamak to visualize.
    """
    fig = plt.figure(figsize=(10,4))
    gs  = gridspec.GridSpec(3, 2)

    gs.update(wspace=0.22, hspace=0, left=0.13, top=0.95, right=0.95, bottom=0.14)

    # Magnetic field
    ds = DREAMSettings()
    tokamak.setMagneticField(ds)

    ax = fig.add_subplot(gs[:,0])
    ds.radialgrid.visualize(ax=ax, ntheta=100)

    # 

    # Density
    r, ne0 = tokamak.getInitialDensity()
    axn = fig.add_subplot(gs[0,1])
    plotSignal(axn, r, ne0/1e19, name=r'$n_e$ ($10^{19}\,$m$^{-3}$)')

    if nticks:
        axn.set_yticks(nticks)

    # Temperature
    r, Te0 = tokamak.getInitialTemperature()
    axT = fig.add_subplot(gs[1,1])
    plotSignal(axT, r, Te0/1e3, name=r'$T_e$ (keV)')

    if Tticks:
        axT.set_yticks(Tticks)

    # Current density
    r, j0 = tokamak.getCurrentDensity()
    axj = fig.add_subplot(gs[2,1])
    plotSignal(axj, r, j0, name=r'$j/j_0$', xlabel=True)

    axj.set_yticks(np.linspace(0.2, 1, 3))

    xticks = np.linspace(0, 1, 6)
    if ds.radialgrid.a == 0:
        xticks *= ds.radialgrid.r_f[-1]
    else:
        xticks *= ds.radialgrid.a

    tickparams = {'direction': 'in', 'length': 5, 'top': True, 'right': True, 'labelsize': 13}
    axn.set_xticks(xticks)
    axn.set_xticklabels([])
    axn.tick_params(**tickparams)
    axT.set_xticks(xticks)
    axT.set_xticklabels([])
    axT.tick_params(**tickparams)
    axj.set_xticks(xticks)
    axj.tick_params(**tickparams)

    def settxt(ax, x, y, txt, direction='tl', **kwargs):
        if 'l' in direction: xdir = 0
        else: xdir = 1
        if 'b' in direction: ydir = 0
        else: ydir = 1

        ax.text(ax.get_xlim()[xdir]*x, ax.get_ylim()[ydir]*y, txt, **kwargs)

    #ax.text(ax.get_xlim()[0]*1.04, ax.get_ylim()[1]*0.92, '(a)')
    #axn.text(
    settxt(ax,  0.86, 0.87, '(a)')
    settxt(axn, 0.88, 0.81, '(b)', 'tr')
    settxt(axT, 0.88, 0.81, '(c)', 'tr')
    settxt(axj, 0.88, 0.81, '(d)', 'tr')

    plt.show()


