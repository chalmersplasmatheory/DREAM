#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMOutput import DREAMOutput

plt.rcParams.update({'font.size': 14})

do_on = DREAMOutput('ware_on_D=0.1.h5')
do_off = DREAMOutput('ware_off_D=0.1.h5')

fig = plt.figure()


#n_hot
n_hot_on = do_on.eqsys.n_hot
n_hot_off = do_off.eqsys.n_hot
plt.plot(do_on.grid.r, do_on.grid.r*n_hot_on[-1,:])
plt.plot(do_off.grid.r, do_off.grid.r*n_hot_off[-1,:],'--')
plt.legend(['ware_on','ware_off'])
plt.xlabel('Radius r')
plt.ylabel('r*n_hot')
plt.title('radius*n_hot')
plt.tight_layout()
plt.show()


#f_hot
fig2 = plt.figure()
do_on.eqsys.f_hot.plot2D(t=-1,r=0,levels=np.linspace(12,24,40))
plt.title('f_hot ware on')
plt.tight_layout()
plt.show()

fig3 = plt.figure()
do_off.eqsys.f_hot.plot2D(t=-1,r=0,levels=np.linspace(12,24,40))
plt.title('f_hot ware off')
plt.tight_layout()
plt.show()


f_hot_diff = do_on.eqsys.f_hot - do_off.eqsys.f_hot
#r=0
fig4 = plt.figure()
f_hot_diff.plot2D(t=-1,r=0,levels=np.linspace(12,24,40))
plt.title('f_hot difference at plasma center')
plt.tight_layout()
plt.show()

#r=-1
fig5 = plt.figure()
f_hot_diff.plot2D(t=-1,r=-1,levels=np.linspace(12,24,40))
plt.title('f_hot difference at plasma edge')
plt.tight_layout()
plt.show()

#n_hot_ratio
fig6 = plt.figure()
n_hot_on = do_on.eqsys.n_hot
n_hot_off = do_off.eqsys.n_hot
n_hot_ratio = n_hot_on/n_hot_off -1
plt.plot(do_on.grid.r, do_on.grid.r*n_hot_ratio[-1,:])
plt.xlabel('radius r')
plt.ylabel('r*n_hot_ratio')
plt.title('radius*n_hot_ratio')
plt.tight_layout()
plt.show()

print(n_hot_ratio[-1,:])
print(sum(n_hot_ratio[-1,:] ))
print(do_on.grid.r)


