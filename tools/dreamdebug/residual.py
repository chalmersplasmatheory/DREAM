# Helper routines for residual vectors output by the DREAM solver.

import h5py
import matplotlib.pyplot as plt
import numpy as np
import warnings

from .petscmat import _mplcursors_frmt1d


#################
# Check for mplcursors
HASMPLCURSORS = False
try:
    import mplcursors
    HASMPLCURSORS = True
except:
    warnings.warn('The recommended package mplcursors was not found.', ImportWarning)


def cmpres(F1, F2, show=True, log=False, eqsys=None):
    """
    Compare two residual vectors.
    """
    dF = np.abs(F2 / F1 - 1)

    if show:
        plotres(dF, show=show, log=log, eqsys=eqsys)


def loadres(filename):
    """
    Load a residual vector from the given MAT file.
    """
    with h5py.File(filename, 'r') as f:
        F = f['F'][:]

    return F


def plotres(res, *args, log=True, show=True, legend=None, eqsys=None):
    """
    Plot the given residual.
    """
    labels = ['Residual 1']


    if log:
        p = lambda x : plt.semilogy(np.abs(x))
    else:
        p = lambda x : plt.plot(x)

    p(res)

    i = 2
    for arg in args:
        p(arg)
        labels.append('Residual {}'.format(i))
        i += 1

    if HASMPLCURSORS:
        cursor = mplcursors.cursor()
        cursor.connect('add', lambda sel : _mplcursors_frmt1d(sel, eqsys=eqsys))

    if legend:
        plt.legend(legend)
    else:
        plt.legend(labels)

    if show:
        plt.show(block=False)


