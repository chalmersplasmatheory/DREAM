#!/usr/bin/env python3
#
# Plot full distribution function (f_hot + f_re)
# ######

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput


do = DREAMOutput('output.h5')

timeindex = -1
xiindex = -2
#plottype = '2d'
#plottype = 'avg'
plottype = 'last'

if plottype == '2d':
    fhot = do.eqsys.f_hot[timeindex,0,:]
    fre  = do.eqsys.f_re[timeindex,0,:]
elif plottype == 'avg':
    fhot = do.eqsys.f_hot.angleAveraged(t=timeindex, r=0)
    fre  = do.eqsys.f_re.angleAveraged(t=timeindex, r=0)
elif plottype == 'last':
    fhot = do.eqsys.f_hot[timeindex,0,xiindex,:]
    fre  = do.eqsys.f_re[timeindex,0,xiindex,:]


if plottype == '2d':
    levels = np.linspace(-66, 25, 20)
    plt.contourf(do.grid.hottail.p, do.grid.hottail.xi, np.log10(np.abs(fhot)), levels=levels)
    plt.contourf(do.grid.runaway.p, do.grid.runaway.xi, np.log10(np.abs(fre)), levels=levels)
    plt.xlim([0, 4])
    plt.ylim([-1,1])
else:
    plt.semilogy(do.grid.hottail.p, fhot)
    plt.semilogy(do.grid.runaway.p, fre)

plt.show()
