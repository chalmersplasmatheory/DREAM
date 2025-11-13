# Functions for evaluating terms of the ion rate equation

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import solve as linsolve
from scipy.constants import e
from scipy.constants import k as k_B


def solve_inv(ni, ne, Te, acd, scd):
    """
    Solve the coronal equilibrium by inverting the linear system
    numerically.
    """
    A = construct_A(ne=ne, Te=Te, acd=acd, scd=scd)
    A[-1,:] = 1
    
    b = np.zeros((A.shape[0],))
    b[-1] = ni

    return linsolve(A, b)


def solve_all(ni, ne, Te, acd, scd):
    """
    Solve the coronal equilibrium.
    """
    Z = acd.Z
    A = np.zeros((Z+1, Z+1))

    I = lambda idx : 0 if idx==Z else scd(Z0=idx, n=ne, T=Te)
    R = lambda idx : 0 if idx==0 else acd(Z0=idx, n=ne, T=Te)

    nout = np.zeros((Z+1,))

    for l in range(Z+1):
        s1, s2 = 0, 0
        if l > 0:
            for j in range(l):
                x = 1
                for k in range(j+1, l+1):
                    x *= R(k)/I(k-1)
                s1 += x

        if l < Z:
            for j in range(l+1, Z+1):
                x = 1
                for k in range(l,j):
                    x *= I(k)/R(k+1)
                s2 += x

        nout[l] = ni / (1+s1+s2)

    return nout


def construct_A(ne, Te, acd, scd):
    """
    Construct the ion rate equation matrix for the specified element at a
    given electron density and temperature.

    :param ne:      Electron density.
    :param Te:      Electron temperature.
    :param acd:     ADAS ACD (recombination) rates.
    :param scd:     ADAS SCD (ionization) rates.
    """
    Z = acd.Z
    A = np.zeros((Z+1, Z+1))

    I = lambda idx : 0 if idx==Z else scd(Z0=idx, n=ne, T=Te)
    R = lambda idx : 0 if idx==0 else acd(Z0=idx, n=ne, T=Te)

    for i in range(Z+1):
        if i > 0:
            A[i,i-1] = I(i-1)*ne

        A[i,i] = -(I(i) + R(i))*ne

        if i < Z:
            A[i,i+1] = R(i+1)*ne

    return A



