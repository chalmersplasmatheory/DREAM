#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMOutput import DREAMOutput

plt.rcParams.update({'font.size': 14})

do = DREAMOutput('output.h5')

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10,4))

# 1.Distribution function evolution
do.eqsys.f_hot.semilogy(t=[0,4,9], ax=axs[0])
axs[0].set_title('Distribution function')

# 2. Runaway rate
rr = do.other.fluid.runawayRate[:,0]
axs[1].plot(do.grid.t[1:], rr, linewidth=2, color='k')
axs[1].set_xlim([0,do.grid.t[-1]])
axs[1].set_ylim([0, 1.2*np.amax(rr)])
axs[1].set_xlabel('Time $t$ (s)')
axs[1].set_ylabel('Runaway rate (s$^{-1}$)')
axs[1].set_title('Runaway rate')

plt.tight_layout()
plt.show()

