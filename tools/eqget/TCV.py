
import h5py
import numpy as np
import sys
import eqhelpers

try:
    from tcvpyget import TCV as TCVMDS
    tcv = TCVMDS()

    AVAILABLE = True
except:
    AVAILABLE = False


def isAvailable():
    """
    Returns ``True`` if this module can be used to fetch equilibrium data
    on this system.
    """
    return AVAILABLE


def getLUKE(shot, time, filename=None):
    """
    Returns magnetic equilibrium data for the given time of the specified
    TCV shot. If ``filename`` is provided, the data is also saved to the
    named LUKE equilibrium data file.

    :param shot: TCV shot to fetch equilibrium data for.
    :param time: Time to fetch equilibrium for.
    :param filename: Name of file to store data in.
    """
    global tcv

    tcv.openShot(shot)
    eq = tcv.equilibrium(time=time)
    tcv.closeall()

    if filename:
        with h5py.File(filename, 'w') as f:
            f.create_group('equil')

            for key in equil.keys():
                f[f'equil/{key}'] = eq[key]

    return eq


def getShaping(shot, time, filename=None, equil=None):
    """
    Fit shaping parameters to the numerical equilibrium in TCV shot 'shot'
    at time 'time'.
    """
    if equil is None:
        equil = getLUKE(shot=shot, time=time)

    return eqhelpers.parametrize_equilibrium(**equil)


