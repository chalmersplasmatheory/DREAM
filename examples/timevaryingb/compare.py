#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import sys
import warnings

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 16})

sys.path.append('../../py')

from DREAM import DREAMOutput
from pathlib import Path

warnings.filterwarnings('ignore', category=RuntimeWarning)


def plotAveragePitch(axs, filename, *args, pIdx=None, plotacctime=False, **kwargs):
    """
    Plot the average pitch angle as a function of energy from the
    given simulation.
    """
    try:
        do = DREAMOutput(filename)
    except: return None, None, None, None

    gxi = do.grid.hottail.xi[:]
    gp  = do.grid.hottail.p[:]
    Vp  = do.grid.hottail.Vprime[:]
    f   = do.eqsys.f_hot
    t   = do.grid.t[:]

    XI = np.matlib.repmat(gxi, gp.size, 1).T

    Vpf = Vp[0,:,:]*f[:,0,:,:]
    dens = np.sum(Vpf, axis=1)
    avxi = np.sum(XI*Vpf, axis=1) / dens
    avTheta = np.arccos(avxi)

    #ax.plot(gp, avxi, *args, label=label, **kwargs)
    if do.settings.radialgrid.dlnB0dt_x is not None:
        dBdt = do.settings.radialgrid.dlnB0dt_x[0] * do.settings.radialgrid.B0[0]
    else:
        dBdt = 0

    label = r'$\dot{{B}} = {:.1f}\,$T/s'.format(dBdt)
    axs[0,0].plot(gp, avTheta[-1,:], *args, label=label, **kwargs)

    # Time evolution
    if pIdx is None:
        pIdx = int(0.75*(gp.size-1))

    axs[0,1].plot(t, avTheta[:,pIdx], label=label, **kwargs)

    if plotacctime:
        tacc = t[np.where(f[:,0,-1,pIdx] > 0.01*f[-1,0,-1,pIdx])][0]
        axs[0,1].plot([tacc, tacc], [-0.5, 0.5], 'k:')
        axs[1,1].plot([tacc, tacc], [-0.5, 0.5], 'k:')

    return t, gp, avTheta, pIdx


fig, axs = plt.subplots(2, 2, figsize=(10,8))

BASEFILE = 'output_dBdt_0.000.h5'
colors = ['tab:blue', 'r', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'yellow']
t, gp, avTheta, pIdx = plotAveragePitch(axs, BASEFILE, color='k', plotacctime=True)

if gp is None:
    raise Exception("A case with dB/dt = 0 must be run for the comparison script to work.")

i = 0
for fname in Path('.').iterdir():
    if fname.suffix != '.h5' or fname.stem[:12] != 'output_dBdt_' or fname.name == BASEFILE:
        continue

    tn, gpn, avThetan, _ = plotAveragePitch(axs, str(fname), color=colors[i], linestyle='--')
    axs[1,0].plot(gp, avThetan[-1,:]-avTheta[-1,:], color=colors[i], linestyle='--')
    axs[1,1].plot(t, avThetan[:,pIdx]-avTheta[:,pIdx], color=colors[i], linestyle='--')

    i += 1

axs[0,0].set_ylabel(r'$\left\langle\theta\right\rangle$')
axs[0,0].set_title('Average pitch angle')
axs[0,0].set_ylim([0, 0.5])
axs[0,0].set_xlim([0, 50])

axs[0,1].set_ylabel(r'$\left\langle\theta\right\rangle$')
axs[0,1].set_title(f'Time evolution ($p={gp[pIdx]:.1f}mc$)')
axs[0,1].set_ylim([0, 0.5])
axs[0,1].set_xlim([0, t[-1]])

for i in range(2):
    axs[1,i].plot([0, 50], [0, 0], 'k:', linewidth=1.0)
    axs[1,i].set_ylabel(r'$\Delta\left\langle\theta\right\rangle$')
    axs[1,i].set_ylim([-0.03, 0.03])

axs[1,0].set_xlabel('Momentum $p/mc$')
axs[1,0].set_xlim([0, 50])

axs[1,1].set_xlabel('Time (s)')
axs[1,1].set_xlim([0, t[-1]])

fig.tight_layout()

for i in range(2):
    box = axs[0,i].get_position()
    axs[0,i].set_position([box.x0, box.y0+box.height*0.10, box.width, box.height*0.90])

for i in range(2):
    box = axs[1,i].get_position()
    axs[1,i].set_position([box.x0, box.y0, box.width, box.height*0.80])

axs[0,0].legend(loc='upper center', bbox_to_anchor=(1.0, -0.05), fancybox=False, edgecolor='k', ncol=2)

plt.show()


