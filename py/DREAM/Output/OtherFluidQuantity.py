# Base class for "other" kinetic (radius + momentum + time) quantities
#

from . FluidQuantity import FluidQuantity


class OtherFluidQuantity(FluidQuantity):
    

    def __init__(self, name, data, description, grid, output):
        """
        Constructor.
        """
        attr = {'description': description}
        super(OtherFluidQuantity, self).__init__(name=name, data=data, grid=grid, attr=attr, output=output)

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
        return '({}) Other fluid quantity of size NT x NR = {} x {}'.format(self.name, self.data.shape[0], self.data.shape[1])


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]


    def _renormalizeTimeIndexForUnknown(self, t):
        """
        Unknowns have one extra initial time point, so map this quantity's time
        indices to unknown time indices by shifting non-negative indices by +1.

        Supports: None, int, list/tuple/ndarray of ints.
        """
        if t is None:
            return slice(1, None, None)

        if isinstance(t, (int, np.integer)):
            return t+1 if t >= 0 else t

        # sequence / fancy indexing: return a list of shifted integers
        if isinstance(t, (list, tuple, np.ndarray, range)):
            return [(int(i)+1) if int(i) >= 0 else int(i) for i in t]

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

        
