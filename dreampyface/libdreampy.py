# Wrapper functions for the C++ functions

from DREAM import DREAMOutput, DREAMSettings, DREAMException
from . Simulation import Simulation
import libdreampy


def get_current_time(ptr): return libdreampy.get_current_time(ptr)
def get_max_time(ptr): return libdreampy.get_max_time(ptr)
def get_others(ptr): return libdreampy.get_others(ptr)
def get_other_data(ptr, name): return libdreampy.get_other_data(ptr=ptr, name=name)
def get_other_info(ptr, name): return libdreampy.get_other_info(ptr=ptr, name=name)
def get_radius_vector(ptr): return libdreampy.get_radius_vector(ptr)
def get_time_vector(ptr): return libdreampy.get_time_vector(ptr)
def get_unknowns(ptr): return libdreampy.get_unknowns(ptr)
def get_unknown_data(ptr, name): return libdreampy.get_unknown_data(ptr=ptr, name=name)
def get_unknown_info(ptr, name): return libdreampy.get_unknown_info(ptr=ptr, name=name)


def register_callback_iteration_finished(func):
    """
    Register a function to be called whenever the non-linear solver
    in DREAM has advanced by one iteration.

    :param func: A callable Python object (e.g. function).
    """
    libdreampy.register_callback_iteration_finished(func)


def register_callback_timestep_finished(func):
    """
    Register a function to be called whenever DREAM has advanced
    by one timestep.

    :param func: A callable Python object (e.g. function).
    """
    libdreampy.register_callback_timestep_finished(func)


def run(ds):
    """
    Run DREAM through the C++ interface.

    :param ds: DREAMSettings object or dict with settings.
    :return: A DREAMOutput object resulting from the simulation.
    """
    if type(ds) == DREAMSettings:
        ds = ds.todict()
    elif type(ds) != dict:
        raise DREAMException("The input object 'ds' must be of type 'DREAMSettings' or 'dict'.")

    d = libdreampy.run(ds)
    return DREAMOutput(d)


def run_simulation(sim):
    """
    Run a previously constructed DREAM simulation.

    :param sim: Simulation object or capsule containing a pointer to the underlying C++ Simulation object.
    """
    if type(sim) == Simulation:
        return sim.run()
    else:
        return libdreampy.run_simulation(sim)


def setup_simulation(ds, returnraw=False):
    """
    Construct a new Simulation object which can be run.

    :param ds:        DREAMSettings object or dictionary containing settings.
    :param returnraw: Return the raw Python capsule used to encapsulate the pointer to the underlying C++ Simulation object.
    """
    if type(ds) == DREAMSettings:
        ds = ds.todict()
    elif type(ds) != dict:
        raise DREAMException("The input object 'ds' must be of type 'DREAMSettings' or 'dict'.")

    #ds['eqsys']['spi']['ZsDrift'] = [int(0)]
    print(ds['eqsys']['spi']['ZsDrift'])
    s = libdreampy.setup_simulation(ds)
    if returnraw:
        return s
    else:
        return Simulation(s)


