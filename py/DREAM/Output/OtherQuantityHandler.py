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

        if other is not None:
            if 'fluid' in other:
                self.fluid = OtherQuantities(other['fluid'], grid, output)
            if 'hottail' in other:
                self.hottail = OtherQuantities(other['hottail'], grid, output, grid.hottail)
            if 'runaway' in other:
                self.runaway = OtherQuantities(other['runaway'], grid, output, grid.runaway)
            if 'scalar' in other:
                self.scalar = OtherQuantities(other['scalar'], grid, output)

    
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
        
        if self.fluid:
            s += "\n   fluid\n{}".format(self.fluid.tostring(padding=6*' '))
        if self.hottail:
            s += "\n   hottail\n{}".format(self.hottail.tostring(padding=6*' '))
        if self.runaway:
            s += "\n   runaway\n{}".format(self.runaway.tostring(padding=6*' '))
        if self.scalar:
            s += "\n   scalar\n{}".format(self.scalar.tostring(padding=6*' '))

        return s


