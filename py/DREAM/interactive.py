"""
This file defines a number of convenience functions which make it
easier to work interactively with the DREAM Python interface.
"""

from .DREAMOutput import DREAMOutput
from .DREAMException import DREAMException

# Declare global variables
_wholist = []

def setup_interactive(do, glob):
    """
    Sets up an interactive session by defining all unknowns as 
    global variables and assigning to them from the given
    DREAM output. 'do' may be either a 'DREAMOutput' object, or
    a string giving the name of the DREAM output to load.

    do:   DREAMOutput object or name of file to load.
    glob: Global environment seen by the caller (this should
          literally be 'globals()')
    """
    global _wholist

    if type(do) == str:
        do = DREAMOutput(do)
    elif not isinstance(do, DREAMOutput):
        raise DREAMException("Unrecognized type of input parameter. Type: {}".format(type(do)))

    _wholist = list(do.eqsys.keys())

    # Declare unknowns
    for uqn in do.eqsys.keys():
        glob[uqn]   = do.eqsys[uqn]

    # Helper functions
    def index2time(i):
        """Converts the specified time index to a time."""
        return do.grid.t[i]

    def time2index(t):
        """Converts the specified time to a time index."""
        return do.time2index(t)

    # Declare other useful stuff
    glob['code'] = do.code
    glob['grid']  = do.grid
    glob['other'] = do.other
    glob['settings'] = do.settings
    glob['solver'] = do.solver
    glob['timings'] = do.timings
    glob['do'] = do
    glob['time2index'] = time2index
    glob['t2i'] = time2index
    glob['index2time'] = index2time
    glob['i2t'] = index2time

    print('Loaded {} unknowns ({})'.format(len(do.eqsys.keys()), do.getFileSize_s()))
    print(do.grid)
    who()

    if len(do.solver.warnings) > 0:
        print(f'\nThe simulation has {len(do.solver.warnings)} warnings.')
        for w in do.solver.warnings:
            print(f'\x1B[1;33m[WARNING]\x1B[0m {w}')


def who():
    """
    Print a list of variables defined in the loaded DREAMOutput object.
    """
    global _wholist
    print('Unknowns:')
    _wholist.sort(key=str.casefold)
    for i in range(0, len(_wholist)):
        if i == 0:
            print('   {}'.format(_wholist[i]), end="")
        else:
            print(', {}'.format(_wholist[i]), end="")

    print("")


