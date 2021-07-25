# Provides access to the UnknownQuantityHandler in the C++ code.


from . import libdreampy


class UnknownQuantityHandler:
    

    def __init__(self, ptr):
        """
        Constructor.

        :param ptr: Pointer to the C++ Simulation object to use for access to the UnknownQuantityHandler.
        """
        self.ptr = ptr


    def getData(self, name):
        """
        Returns the most recent solution for the named unknown quantity.
        The method returns a dict with the following keys:
          
          'p':  Momentum grid (if applicable)
          'xi': Pitch grid (if applicable)
          'r':  Radial grid (if applicable)
          't':  Time grid
          'x':  Unknown data
        """
        return libdreampy.get_unknown_data(self.ptr, name)[name]


    def getInfo(self, name=None):
        """
        Returns information about one or all unknown quantities in
        the simulation.

        :param name: Name of unknown quantity to return info for. If ``None``, return info for all quantities.
        """
        if name is None:
            return libdreampy.get_unknowns(self.ptr)
        else:
            return libdreampy.get_unknown_info(self.ptr)

