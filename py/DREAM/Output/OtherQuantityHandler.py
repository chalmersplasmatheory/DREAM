# Thin wrapper for other quantities

from . OtherQuantities import OtherQuantities


class OtherQuantityHandler:
    

    def __init__(self, other=None, grid=None, output=None):
        """
        Constructor.
        """
        self.fluid   = None
        self.hottail = None
        self.runaway = None
        self.scalar  = None

        self.categories = []

        if other is not None:
            for category in other.keys():
                self.categories.append(category)
                setattr(self, category, OtherQuantities(category, other[category], grid, output))

    
    def __contains__(self, item):
        return (item in self.__dict__)


    def __getitem__(self, index):
        return self.__dict__[index]


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        s = "OtherQuantityHandler with"
        
        for category in self.categories:
            s += "\n   {}\n{}".format(category, self[category].tostring(padding=6*' '))

        return s


