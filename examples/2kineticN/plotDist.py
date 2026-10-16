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

plt.figure()

# Plot electron density evolution
nhot = do.eqsys.f_hot.density(r=0)
nre  = do.eqsys.f_re.density(r=0)
n0   = do.eqsys.n_cold[0,0]
plt.plot(do.grid.t, (nhot+nre-n0)/n0, 'k', linewidth=3)
plt.ylabel(r'$(n_{\rm hot} + n_{\rm RE}) / n_{\rm e,0}$')
#plt.semilogy(do.grid.t, nhot, '--', linewidth=2)
#plt.semilogy(do.grid.t, nre, '--', linewidth=2)
#plt.legend([r'$n_{\rm e}$', r'$n_{\rm hot}$', r'$n_{\rm RE}$'])

plt.show()
