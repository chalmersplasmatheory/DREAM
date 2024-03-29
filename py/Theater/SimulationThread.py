# A class for running a DREAM simulation on a separate thread

from PyQt5.QtCore import QThread, pyqtSignal
from dreampyface import Simulation
from DREAM import DREAMOutput


class SimulationThread(QThread):
    

    timestepTaken = pyqtSignal(Simulation)
    output = None


    def __init__(self, simulation):
        """
        Constructor.
        """
        super().__init__()

        self.simulation = simulation


    def getOutput(self):
        """
        Returns the DREAMOutput object resulting from the simulation.
        """
        return self.output


    def run(self):
        """
        Run the simulation.
        """
        self.simulation.onTimestep(lambda sim : self.timestepCompleted(sim))
        self.output = self.simulation.run()


    def timestepCompleted(self, sim):
        """
        Function to be called whenever the simulation has completed
        another time step.
        """
        self.timestepTaken.emit(sim)


