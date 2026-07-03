#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMOutput import DREAMOutput

plt.rcParams.update({'font.size': 14})

do_on = DREAMOutput('t=4e-2/ware_on.h5')
do_off = DREAMOutput('t=4e-2/ware_off.h5')

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,10))


#n_hot
n_hot_on = do_on.eqsys.n_hot
n_hot_off = do_off.eqsys.n_hot
axs[0][0].plot(do_on.grid.r, do_on.grid.r*n_hot_on[-1,:])
axs[0][0].plot(do_off.grid.r, do_off.grid.r*n_hot_off[-1,:])
axs[0][0].legend(['ware_on','ware_off'])
axs[0][0].set_xlabel('Radius r')
axs[0][0].set_ylabel('r*n_hot')
axs[0][0].set_title('radius*n_hot')


#n_re
#n_re_on = do_on.eqsys.n_re
#n_re_off = do_off.eqsys.n_re
#axs[0][1].plot(do_on.grid.r, do_on.grid.r*n_re_on[-1,:])
#axs[0][1].plot(do_off.grid.r, do_off.grid.r*n_re_off[-1,:])
#axs[0][1].legend(['ware on','ware off'])
#axs[0][1].set_xlabel('Radius r')
#axs[0][1].set_ylabel('r*n_re')
#axs[0][1].set_title('radius*n_re')

#f_hot
do_on.eqsys.f_hot.plot2D(t=-1,r=0,levels=np.linspace(12,24,40), ax=axs[1][0])
axs[1][0].set_title('f_hot ware on')

do_off.eqsys.f_hot.plot2D(t=-1,r=0,levels=np.linspace(12,24,40), ax=axs[1][1])
axs[1][1].set_title('f_hot ware off')
#plt.colorbar(ax =axs[1][1])


f_hot_diff = do_on.eqsys.f_hot - do_off.eqsys.f_hot
f_hot_diff.plot2D(t=-1,r=0,levels=np.linspace(12,24,40), ax=axs[1][2])
axs[1][2].set_title('f_hot ware_on-ware_off')

f_hot_diff.plot2D(t=-1,r=-1,levels=np.linspace(12,24,40), ax=axs[0][1])
axs[0][1].set_title('f_hot ware_on-ware_off')

#n_hot_ratio
n_hot_on = do_on.eqsys.n_hot
n_hot_off = do_off.eqsys.n_hot
n_hot_ratio = n_hot_on/n_hot_off -1
axs[0][2].plot(do_on.grid.r, do_on.grid.r*n_hot_ratio[-1,:])
axs[0][2].set_xlabel('Radius r')
axs[0][2].set_ylabel('r*n_hot_ratio')
axs[0][2].set_title('radius*n_hot_ratio')

print(len(do_on.grid.t))
print(do_on.grid.t[-1])

plt.tight_layout()
plt.show()


