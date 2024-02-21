# Provides access to the OtherQuantityHandler in the C++ code.


from . import libdreampy


class OtherQuantityHandler:
    

    def __init__(self, ptr):
        """
        Constructor.

        :param ptr: Pointer to the C++ Simulation object to use for access to the OtherQuantityHandler.
        """
        self.ptr = ptr


    def getData(self, name):
        """
        Returns the most recent solution for the named other quantity.
        The method returns a dict with the following keys:

          'p':  Momentum grid (if applicable)
          'xi': Pitch grid (if applicable)
          'r':  Radial grid (if applicable)
          't':  Time grid
          'x':  Other quantity data
        """
        return libdreampy.get_other_data(self.ptr, name)[name]


    def getInfo(self, name=None):
        """
        Returns information about one or all other quantities in
        the simulation. Each entry contains the following fields:

          description: Description of unknown quantity.
          nelements:   Total number of elements used for other quantity.
          nmultiples:  Number of "multiples" represented by the other quantity (e.g. number of ion species and charge states).
        """
        if name is None:
            return libdreampy.get_others(self.ptr)
        else:
            return libdreampy.get_other_info(self.ptr, name)

