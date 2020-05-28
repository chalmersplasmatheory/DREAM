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

        if other is not None:
            if 'fluid' in other:
                self.fluid = OtherQuantities(other['fluid'], grid, output)
            if 'hottail' in other:
                self.hottail = OtherQuantities(other['hottail'], grid, output, grid.hottail)
            if 'runaway' in other:
                self.runaway = OtherQuantities(other['runaway'], grid, output, grid.runaway)

    
    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        l = []

        if self.fluid is not None: l.append('fluid')
        if self.hottail is not None: l.append('hottail')
        if self.runaway is not None: l.append('runaway')

        if len(l) > 1:
            s = ', '.join(l[:-1]) + " and " + l[-1]
        else:
            s = l[0]

        return "OtherQuantityHandler with {}".format(s)


