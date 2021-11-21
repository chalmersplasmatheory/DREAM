# Implementation of a data provider encapsuling a dreampyface Simulation
# object.


from . DataProvider import DataProvider


class DataProviderSimulation(DataProvider):
    

    def __init__(self, simulation):
        """
        Constructor.

        :param simulation: Simulation object to encapsulate.
        """
        self.simulation = simulation


    def getRadius(self):
        """
        Returns the radius vector of the data.
        """
        return self.simulation.getRadius()


    def getTime(self):
        """
        Returns the time vector of the data.
        """
        return self.simulation.getTime()


    def getOtherInfo(self, name=None):
        """
        Returns basic information about the other quantities
        of the equation system.
        """
        return self.simulation.others.getInfo(name=name)


    def getUnknownInfo(self, name=None):
        """
        Returns basic information about the unknowns of the
        equation system.
        """
        return self.simulation.unknowns.getInfo(name=name)
    
