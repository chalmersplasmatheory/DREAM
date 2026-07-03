#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMOutput import DREAMOutput

plt.rcParams.update({'font.size': 14})

do = DREAMOutput('t=5e-3_Nt=100/ware_on.h5')

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10,4))

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


Ar = do.other.hottail.Ar
print(Ar)
print(Ar[-1,:,33,0])
axs[2].plot(do.grid.t[1:], Ar[:,3,33,53])
axs[2].set_xlabel('Time $t$ (s)')
axs[2].set_ylabel('Ar')
axs[2].set_title('Advection Coefficient')


#print(np.sum(Ar[:,5,33,53]))


plt.tight_layout()
plt.show()

