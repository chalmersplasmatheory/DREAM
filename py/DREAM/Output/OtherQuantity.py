# Base class for an "other" quantity which is neither kinetic, fluid nor scalar.

from . UnknownQuantity import UnknownQuantity


class OtherQuantity(UnknownQuantity):
    

    def __init__(self, name, data, description, grid, output):
        """
        Constructor.
        """
        attr = {'description': description}
        super().__init__(name=name, data=data, grid=grid, attr=attr, output=output)

        self.time = grid.t[1:]


    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        #s = self.__str__() 
        return self.__str__()


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Other quantity of size {}'.format(self.name, self.data.shape)


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def _renormalizeTimeIndexForUnknown(self, t):
        """
        Tries to re-normalize the given time index so that it correctly indexes
        a regular unknown quantity (which has a different time base).
        """
        if t is None:
            t = slice(None)

        start, stop, step = t.start, t.stop, t.step
        if start is None or start == 0:
            start = 1
        if stop is not None:
            stop += 1

        t = slice(start, stop, step)

        return t

    def new_like(self, name=None, data=None, grid=None, output=None, attr=None, description=None):
        """
        Creates a new object of the same type where the provided quantities replace
        those of self.
        """
        if name is None:
            name = self.name
        if data is None:
            data = self.data
        if grid is None:
            grid = self.grid
        if output is None:
            output = self.output
        if description is None:
            if attr is not None and "description" in attr:
                description = attr["description"]
            else:
                description = self.description
        
        return type(self)(name=name, data=data, grid=grid, output=output, description=description)

    def getMultiples(self):
        """
        Get the number of "multiples" (e.g. number of ion species and
        charge states) covered by this quantity. The total number of elements
        in 'self.data' is the size of the grid on which this quantity lives
        (i.e. scalar grid, fluid grid, or a kinetic grid) times this number.
        """
        return 1

        
