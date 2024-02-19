# Various high-level functions for running DREAM

from . import resolvedreampaths
from PyQt5 import QtWidgets
import numpy as np
import sys
import time

from DREAM import DREAMSettings
from dreampyface import Simulation

from . PlotWindow import PlotWindow


def plot_fluidquantity(simulation, unknown, radius=0):
    """
    Plot the fluid quantity named ``unknown`` at the given radius.
    The parameter ``radius`` may be an array, in which case each
    radius specified is plotted. If ``radius`` is ``None``, all radii
    are plotted.

    :param simulation: dreampyface.Simulation object to fetch data from.
    :param unknown:    Name of unknown quantity to plot.
    :param radius:     Indices for radial points to plot.
    """
    u = simulation.unknowns.getData(unknown)

    # Select all radii?
    if radius is None:
        radius = u['r']

    # Plot just one radius?
    if np.isscalar(radius):
        radius = [radius]

    x, y = [], []
    for r in radius:
        x.append(u['t'])
        y.append(u['x'][:,r])
    
    return x, y


def run_monitor(s, func, *args, **kwargs):
    """
    Run a DREAM simulation and monitor the evolution of the named
    unknown while doing so.

    :param s:      DREAMSettings, dict or dreampyface.Simulation object.
    :param func:   Function returning data to monitor. Should take a
                   dreampyface.Simulation object as input and return two lists
                   containing the x and y vectors to plot (i.e. each list contains
                   one or more vectors of data).
    :param args:   Passed on to 'PlotWindow'.
    :param kwargs: Passed on to 'PlotWindow'.
    """
    if type(s) == dict:
        s = DREAMSettings(s)
    if type(s) == DREAMSettings:
        s = Simulation(s)
    if type(s) != Simulation:
        raise Exception("The parameter 's' must be either a settings dict, a DREAM.DREAMSettings object or a dreampyface.Simulation object.")

    # TODO run
    app = QtWidgets.QApplication(sys.argv)
    plotWindow = PlotWindow(simulation=s, callback=func, *args, **kwargs)

    plotWindow.show()
    plotWindow.run()
    app.exec_()

    return plotWindow.getOutput()


