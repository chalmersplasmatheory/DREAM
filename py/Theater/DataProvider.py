# Base class defining the interface to classes which provide data for the
# Theater GUI.


class DataProvider:
    

    def __init__(self):
        """
        Constructor.
        """
        pass


    def getRadius(self):
        """
        Returns the radius vector of the data.
        """
        raise NotImplementedError()


    def getTime(self):
        """
        Returns the time vector of the data.
        """
        raise NotImplementedError()


    def getOtherInfo(self, name=None):
        """
        Get basic information about all other quantities.
        """
        raise NotImplementedError()


    def getUnknownInfo(self, name=None):
        """
        Get basic information about all unknowns.
        """
        raise NotImplementedError()


