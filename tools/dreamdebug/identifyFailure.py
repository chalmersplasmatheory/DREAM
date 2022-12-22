#!/usr/bin/env python3
#
# IDENTIFY UNKNOWN QUANTITY CAUSING DIVERGENCE
#
# Sometimes DREAM simulations fail to converge due to errors in the jacobian
# matrix. Usually, the solver ends up iterating with multiple quantities
# failing to converge. This script is intended to be able to identify which of
# the various unknown quantities in a system of equations that is actually
# causing the divergence.
#
# The script takes a DREAM jacobian matrix as input and iterates through the
# various unknowns of the equation system, fixing the value of the unknown
# and attempting to solve the resulting system of equations. In principle,
# the unknown which, when given a fixed value, yields the smallest variation
# in the solution should be the quantity causing the problems.
#
# HOW TO USE
#
# 1. Run DREAM and export a (i) jacobian, (ii) residual and (iii) solution
#    for a time step/iteration which seems problematic.
#
# 2. Run this script, specifying in which directory the files in step 1 are
#    located. Also, if more than one timestep and iteration is availalbe, you
#    can specify which time step/iteration to use using the '-t' and '-i'
#    flags, respectively.
#
#    ./IdentifyFailure.py debugdir -t 1 -i 2
#

import argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import sys

if 'DREAMPATH' in os.environ:
    sys.path.append(os.environ['DREAMPATH'])

from DREAM import DREAMOutput
from scipy.sparse import csr_matrix, find
from tools.dreamdebug import DREAMEqsys, petscmat
from pathlib import Path


ITERATION = None
TIMESTEP = None


def getFilename(basename, debugdir, timestep, iteration, suffix='', returnTimestep=False):
    """
    Determine the filename of a file with general format

      {debugdir}/{basename}_{timestep}_{iteration}.{suffix}
    """
    global TIMESTEP, ITERATION

    s = ''
    if debugdir:
        s += debugdir + '/'

    # Automatically determine timestep (select first timestep)
    if timestep is None or iteration is None:
        if TIMESTEP is not None and ITERATION is not None:
            timestep = TIMESTEP
            iteration = ITERATION
        else:
            b = f'{basename}_'
            if timestep is None: b += '[0-9]*_'
            else: b += f'{timestep}_'

            if iteration is None: b += '[0-9*]'
            else: b += f'{iteration}'

            ptrn = re.compile(b)
            ts = []
            it = []
            for p in Path(debugdir).iterdir():
                m = ptrn.match(p.stem)
                if m:
                    ss = p.stem.split('_')
                    ts.append(int(ss[-2]))
                    it.append(int(ss[-1]))

            if timestep is None:
                TIMESTEP = timestep = min(ts)
            if iteration is None:
                ITERATION = iteration = min(it)

            print(f'Selecting timestep {timestep}, iteration {iteration}.')

    s += f'{basename}_{timestep}_{iteration}'

    if suffix:
        if suffix[0] == '.':
            s += f'{suffix}'
        else:
            s += f'.{suffix}'

    p = Path(s).resolve().absolute()

    if returnTimestep:
        return str(p), timestep, iteration
    else:
        return str(p)


def loadEqsys(nions, filename='eqsys.txt'):
    """
    Try to load the equation system specification.
    """
    try:
        return DREAMEqsys(filename, nions=nions)
    except:
        return None


def loadJacobianAndResidual(debugdir, timestep=None, iteration=None):
    """
    Load the jacobian matrix for the specified time step and iteration.
    """
    global ITERATION, TIMESTEP
    
    if timestep is None:
        timestep = TIMESTEP
    if iteration is None:
        iteration = ITERATION

    Aname = getFilename('petsc_jac', debugdir, timestep, iteration)
    bname = getFilename('residual', debugdir, timestep, iteration, suffix='mat')
    A = petscmat.load(Aname)
    with h5py.File(bname) as f:
        b = f['F'][:]

    return A, b


def loadOutput(debugdir, timestep, iteration):
    """
    Load a DREAM debug output file.
    """
    fname = getFilename('debugout', debugdir, timestep, iteration, suffix='h5')
    return DREAMOutput(fname, loadsettings=False)


def calculateVariation(eqsys, do, A, b):
    """
    Calculate the variation in the solution when different unknowns are
    given fixed values.
    """
    i, x, y = 0, 0, 0
    n = len(eqsys.unknowns)
    arr = np.zeros((n, n))
    namelist = []

    for u in eqsys.unknowns:
        namelist.append(u.name)

        # Construct new equation system
        A2 = A.copy()
        b2 = b.copy()

        # Set equation corresponding to current unknown to
        #   x = constant
        A2[i:i+u.size,:] = 0
        b2[i:i+u.size] = do.eqsys[u.name].data[:].flatten()

        # (add ones along block diagonal for identity term)
        for j in range(u.size):
            A2[i+j,i+j] = 1

        i += u.size

        # Solve modified equation system
        solution = petscmat.solve(A2, b2)

        # Evaluate how much all other unknowns have changed
        # in the modified solution
        k, y = 0, 0
        for u2 in eqsys.unknowns:
            if np.sum(do.eqsys[u2.name].data[:]) == 0:
                arr[y,x] = 0
            else:
                arr[y,x] = np.abs(np.sum(solution[k:k+u2.size]) / np.sum(do.eqsys[u2.name].data[:]))

            k += 1
            y += 1

        x+= 1

    # Estimate induced variation by summing over all columns in each row.
    # (we do this row by row to handle unknowns which are 0 in a special way)
    colsum = np.zeros(n)
    for col in range(n):
        if np.linalg.norm(do.eqsys[eqsys.unknowns[col].name].data[:]) == 0:
            colsum[col] = np.inf
        else:
            colsum[col] = np.sum(arr[:,col])

    # Unknown with smallest variation is likely the troublesome unknown
    return namelist, colsum


def printVariation(namelist, variation):
    """
    Print the result of the calculation.
    """
    i = 1
    for name, val in sorted(zip(namelist, variation), key=lambda x : x[1]):
        print(f'{i:02d}. {name:12s} = {val:.12e}')

        i += 1


def parse_args():
    parser = argparse.ArgumentParser("DREAM debug")

    parser.add_argument('-i', '--iteration', help="Iteration number to load", type=int, nargs='?', default=None)
    parser.add_argument('-t', '--timestep', help="Time step to load", type=int, nargs='?', default=None)
    parser.add_argument('debugdir', help="Directory containing debug output", nargs='?', default='.')

    return parser.parse_args()


def main():
    """
    Program entry point.
    """
    global ITERATION, TIMESTEP

    args = parse_args()

    do = loadOutput(args.debugdir, timestep=args.timestep, iteration=args.iteration)
    eqsys = loadEqsys(nions=len(do.eqsys.n_i.ions))

    A, b = loadJacobianAndResidual(args.debugdir)

    namelist, variation = calculateVariation(eqsys, do, A, b)

    printVariation(namelist, variation)

    return 0


if __name__ == '__main__':
    sys.exit(main())


