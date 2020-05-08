"""
This file defines a number of convenience functions which make it
easier to work interactively with the DREAM Python interface.
"""

from .DREAMOutput import DREAMOutput
from .DREAMException import DREAMException

# Declare global variables
_wholist = []
E_field = None
n_cold  = None
n_hot   = None
n_re    = None
T_cold  = None
n_tot   = None
n_i     = None
f_hot   = None

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
    global _wholist, E_field, f_hot, n_cold, n_hot, n_re, n_tot, n_i, T_cold

    if type(do) == str:
        do = DREAMOutput(do)
    elif type(do) != DREAMOutput:
        raise DREAMException("Unrecognized type of input parameter.")

    _wholist = list(do.eqsys.keys())

    # Declare unknowns
    glob['E_field'] = do.eqsys.E_field
    glob['n_cold']  = do.eqsys.n_cold
    glob['n_hot']   = do.eqsys.n_hot
    glob['n_re']    = do.eqsys.n_re
    glob['T_cold']  = do.eqsys.T_cold
    glob['n_tot']   = do.eqsys.n_tot
    glob['n_i']     = do.eqsys.n_i
    glob['f_hot']   = do.eqsys.f_hot

    print('Loaded {} unknowns.'.format(len(do.eqsys.keys())))
    print(do.grid)
    who()


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


