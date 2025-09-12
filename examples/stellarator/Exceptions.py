class BadSimulationException(Exception):
    """ This DREAM simulation is bad. """
    pass

class MaximumTimestepsException(BadSimulationException):
    """ Maximum DREAM time resolution was reached. """
    pass

class MinimumStepLengthException(BadSimulationException):
    """ Minimum DREAM time step length was reached. """
    pass

class TransportException(BadSimulationException):
    """ No thermal quench with this initial heat transport. """
    pass

class DiscardSimulation(BadSimulationException):
    """ Simulation crashed. """
    pass

class timeOutException(BadSimulationException):
    """ Simulation crashed. """
    pass

class QuitException(Exception):
    """ This DREAM simulation is bad. """
    pass