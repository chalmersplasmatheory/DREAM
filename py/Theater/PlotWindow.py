
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtCore import QThread, pyqtSignal

import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import numpy as np
from dreampyface import Simulation

from DREAM import GeriMap


mpl.rc('text', usetex=False)


class PlotWindow(QtWidgets.QFrame):


    def __init__(self, simulation, callback, xlim=(0, 1.0), ylim=(0, 1.0),
                 width=600, height=600, parent=None):
        """
        Constructor.

        :param simulation: Simulation settings or dreampyface.Simulation object to run.
        :param callback:   Function which retrieves data to plot. Takes a dreampyface.Simulation object as input and returns a time vector and a data vector.
        """
        super().__init__(parent)

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.setWindowTitle('Plot window')

        if type(simulation) != Simulation:
            raise Exception("Expected 'Simulation' to be of type dreampyface.Simulation.")

        self.simthread = None
        self.simulation = simulation
        self.callback = callback
        self.output = None

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.resize(width,height)

        self.ax = self.figure.add_subplot(111)
        self.lines = None

        self.xlim = xlim
        self.ylim = ylim

        x, y = self.callback(self.simulation)
        self.updatePlot(x, y)
        plt.pause(0.05)


    def drawSafe(self):
        try:
            self.canvas.draw()
        except RuntimeError as e:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText(e.strerror)
            msg.setWindowTitle('Runtime Error')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()


    def getOutput(self):
        """
        Returns the 'DREAMOutput' object produced by the simulation.
        If the Simulation has not finished, returns 'None'.
        """
        return self.output


    def run(self):
        """
        Run the simulation.
        """
        if self.simthread != None:
            raise Exception("Cannot run more than one simulation at a time.")

        self.simthread = SimulationThread(self.simulation, self.callback)
        self.simthread.timestepTaken.connect(self.updatePlot)
        self.simthread.start()


    def updatePlot(self, x, y):
        """
        Update the current plot.
        """
        N = len(x)
        colors = GeriMap.get(N=N+1)

        if self.lines is None:
            self.lines = []
            for i in range(N):
                clr = colors(i/(N+1))
                l, = self.ax.plot(x[i], y[i], color=clr)
                self.lines.append(l)

            self.ax.set_xlim(self.xlim)
            self.ax.set_ylim(self.ylim)
        else:
            for i in range(N):
                l = self.lines[i]
                l.set_data(x[i], y[i])

        self.canvas.draw()
        self.canvas.flush_events()


class SimulationThread(QThread):
    

    timestepTaken = pyqtSignal(list, list)


    def __init__(self, simulation, callback):
        """
        Constructor.
        """
        super().__init__()

        self.callback = callback
        self.simulation = simulation


    def run(self):
        """
        Run the simulation.
        """
        self.simulation.onTimestep(lambda ptr : self.timestepCompleted(ptr))
        out = self.simulation.run()


    def timestepCompleted(self, ptr):
        """
        Function to be called whenever the simulation has completed
        another time step.
        """
        s = Simulation(ptr)
        x, y = self.callback(s)

        self.timestepTaken.emit(x, y)


