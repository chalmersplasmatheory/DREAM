# Thin wrapper for static quantities

from .StaticQuantities import StaticQuantities


class StaticQuantityHandler:
    

    def __init__(self, static=None):
        """
        Constructor.
        """
        self.categories = []

        if static is not None:
            for category in static.keys():
                self.categories.append(category)
                setattr(self, category, StaticQuantities(category, static[category]))


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
        s = "StaticQuantityHandler with"

        for category in self.categories:
            s += f"\n   {category}\n{self[category].tostring(padding=6*' ')}"

        return s


