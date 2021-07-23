# Definition of the pyface 'Simulation' object which can be used
# to interact with a C++ 'Simulation' object, used internally by DREAM.

from . import libdreampy
from DREAM import DREAMSettings, DREAMOutput
from . UnknownQuantityHandler import UnknownQuantityHandler


class Simulation:
    

    def __init__(self, ptr):
        """
        Constructor.

        :param ptr: Capsule object containing a pointer to the C++ 'Simulation' object.
        """
        if type(ptr) == str:
            self.ptr = libdreampy.setup_simulation(DREAMSettings(ptr), returnraw=True)
        elif type(ptr) == DREAMSettings:
            self.ptr = libdreampy.setup_simulation(ptr, returnraw=True)
        else:
            # ...otherwise we assume it's a Python capsule, containing
            # a pointer to the underlying C++ Simulation object
            self.ptr = ptr


        self.unknowns = UnknownQuantityHandler(self.ptr)
    

    def getCurrentTime(self):
        """
        Returns the current time of the simulation.
        """
        return libdreampy.get_current_time(self.ptr)


    def getMaxTime(self):
        """
        Returns the maximum simulation time.
        """
        return libdreampy.get_max_time(self.ptr)


    def run(self):
        """
        Run this simulation.
        """
        return DREAMOutput(libdreampy.run_simulation(self.ptr))


