# Base class for a static quantity


class StaticQuantity:
    

    def __init__(self, name, data, description=''):
        """
        Constructor.
        """
        self.name = name
        self.data = data
        self.description = description


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return f"({self.name}) Static quantity of size {self.data.shape}"


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


