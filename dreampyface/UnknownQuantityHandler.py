# Provides access to the UnknownQuantityHandler in the C++ code.


from . import libdreampy


class UnknownQuantityHandler:
    

    def __init__(self, ptr):
        """
        Constructor.

        :param ptr: Pointer to the C++ Simulation object to use for access to the UnknownQuantityHandler.
        """
        self.ptr = ptr


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

