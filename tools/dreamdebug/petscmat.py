#
# Routines for loading and viewing sparse PETSc matrices stored in
# the 'BINARY_MATLAB' format. These routines also provide support for
# DREAMEqsys objects, which allows more informative markers to be placed.

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.sparse as sparse
import scipy.sparse.linalg
import sys
import warnings

#################
# Check for mplcursors
HASMPLCURSORS = False
try:
    import mplcursors
    HASMPLCURSORS = True
except:
    warnings.warn('The recommended package mplcursors was not found.', ImportWarning)

matplotlib.rcParams.update({'text.usetex': False})
warnings.filterwarnings('ignore', category=sparse.SparseEfficiencyWarning)

#################
# Load PETSc routine for reading 'BINARY_MATLAB' files
if 'PETSC_DIR' not in os.environ:
    raise Exception("Unable to locate PETSc. Please set the 'PETSC_DIR' environment variable.")
else:
    sys.path.append(os.environ['PETSC_DIR'] + '/lib/petsc/bin')

import PetscBinaryIO


def cmp(m1, m2, show=True, tollow=None, tolup=None, eqsys=None):
    """
    Compare two matrices

    m1, m2: Matrices to compare.
    tollow: Lower tolerance for an element to be considered different.
    tolup:  Upper tolerance for an element to be considered "extremely" different.
    """
    r, c = sparse.find(m1)[:-1]

    dm  = m1.copy()
    dm[r,c] = np.abs(m2[r,c] / m1[r,c] - 1)

    if tollow is not None:
        r, c, _ = sparse.find(dm > tollow)
        dm2 = sparse.csr_matrix(dm.shape)
        dm2[r,c] = dm[r,c]
        dm = dm2

    spy(dm, show=False, markersize=1, eqsys=eqsys)

    if tolup is not None:
        r, c, _ = sparse.find(dm > tolup)
        plt.plot(c, r, 'ro')

    if show:
        plt.show(block=False)

    return dm


def load(filename):
    """
    Loads the named PETSc binary matrix.
    """
    pbio = PetscBinaryIO.PetscBinaryIO()
    mat = pbio.readBinaryFile(filename)[0]
    return sparse.csr_matrix((mat[1][2], mat[1][1], mat[1][0]))


def _mplcursors_frmt1d(sel, eqsys=None):
    """
    Matplotlib cursor format for regular 1D plots.
    """
    x, y = sel.target

    rx = int(np.round(x))

    sx = 'x: {:d}'.format(rx)
    sy = 'y: {:.5e}'.format(y)

    if eqsys is not None:
        sx += ' ({})'.format(eqsys[rx].toxstring(rx))

    s = sx + '\n' + sy

    sel.annotation.set_text(s)


def _mplcursors_frmt2d(sel, mat=None, eqsys=None):
    """
    Matplotlib cursor format for spy() on sparse matrices.
    """
    c, r = sel.target
    
    sr = 'Row: {:d}'.format(int(r))
    sc = 'Col: {:d}'.format(int(c))

    if eqsys is not None:
        sr += ' ({})'.format(eqsys[r].toxstring(r))
        sc += ' ({})'.format(eqsys[c].toxstring(c))

    s = sr + '\n' + sc

    if mat is not None:
        sv = 'Val: {:.4e}'.format(mat[r,c])
        s += '\n' + sv

    sel.annotation.set_text(s)


def plotcol(m1, m2=None, col=None, log=False, show=True, legend=None, eqsys=None):
    """
    Plot the specified rows of the matrices.
    """
    if col is None:
        raise Exception("Parameter 'col' not specified.")

    if log:
        plt.semilogy(np.abs(m1[:,col].A))
        if m2 is not None: plt.semilogy(np.abs(m2[:,col].A))
    else:
        plt.plot(m1[:,col].A)
        if m2 is not None: plt.plot(m2[:,col].A)

    if HASMPLCURSORS:
        cursor = mplcursors.cursor()
        cursor.connect('add', lambda sel : _mplcursors_frmt1d(sel, eqsys=eqsys))

    if legend:
        plt.legend(legend)
    else:
        if m2 is not None: plt.legend(['Matrix 1', 'Matrix 2'])

    if show:
        plt.show(block=False)


def plotcoll(*args, **kwargs):
    plotcol(*args, log=True, **kwargs)


def plotrow(m1, m2=None, row=None, *args, **kwargs):
    if m2 is None:
        plotcol(m1.T, None, col=row, *args, **kwargs)
    else:
        plotcol(m1.T, m2.T, col=row, *args, **kwargs)


def plotrowl(*args, **kwargs):
    plotrow(*args, log=True, **kwargs)


def setxlabels(dreameqsys, ax=None):
    """
    Set x labels of current plot based on the given DREAMEqsys object.
    """
    if ax is None:
        ax = plt.gca()

    ax.set_xticks(dreameqsys.getoffsets())
    ax.set_xticklabels(dreameqsys.getnames())

    plt.gcf().canvas.draw()


def solve(A, b, **kwargs):
    """
    Solves the sparse linear system

      Ax = b
    
    and returns the vector x.
    """
    return scipy.sparse.linalg.spsolve(A, b, **kwargs)


def spy(m, show=True, aspect='auto', eqsys=None, *args, **kwargs):
    """
    Spy on sparse matrix.
    """
    global HASMPLCURSORS
    x = plt.spy(m, aspect=aspect, *args, **kwargs)

    if HASMPLCURSORS:
        cursor = mplcursors.cursor()
        cursor.connect('add', lambda sel : _mplcursors_frmt2d(sel, mat=m, eqsys=eqsys))

    if show:
        plt.show(block=False)

    return x


