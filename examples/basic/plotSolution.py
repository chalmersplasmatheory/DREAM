#!/usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import numpy as np


E_field = None
n_hot = None
n_re = None
f_hot = None
n_cold = None
n_i = None
with h5py.File('output.h5', 'r') as f:
    E_field = f['E_field/x'][:]
    n_hot = f['n_hot/x'][:]
    n_re = f['n_re/x'][:]
    f_hot = f['f_hot/x'][:]
    n_cold = f['n_cold/x'][:]
    n_i = f['n_i/x'][:]


def plotFhot(it=0, ir=0, ixi=0):
    global f_hot

    plt.plot(f_hot[it,ir,ixi,:])
    plt.show()


def plotNhot(it=0):
    global n_hot

    plt.plot(n_hot[it,:])
    plt.show()

def plotNion(it=0):
    global n_i

    plt.plot(n_i[it,:])
    plt.show()

#plotFhot(1)
#plotNhot(0)
plotNion(1)

