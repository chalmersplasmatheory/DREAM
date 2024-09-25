#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import e, atomic_mass

plt.rcParams.update({
    'font.family': 'FreeSans',
    'font.size': 15
})

from DREAM import DREAMOutput


do = DREAMOutput('test.h5')

TIMESTEP = 0
gamma = 7
n = do.eqsys.n_cold[:,-1]
T = do.eqsys.T_cold[:,-1]
mi = 3.3435837724e-27
cs = np.sqrt(e * T / mi)
qPar = gamma * n * cs * e*T

print(f'T = {e*T[TIMESTEP]}')
print(f'mi = {mi}')
print(f'Fr = {gamma * n[TIMESTEP] * cs[TIMESTEP] * e}')
#print(qPar)

plt.semilogy(do.grid.t[:-1], do.other.scalar.Wcold_Tcold_Ar[:,0])
plt.semilogy(do.grid.t, gamma*n*cs*e)

plt.figure(figsize=(6,4))

plt.semilogy(do.grid.t, qPar, lw=3, label=r'$q_\parallel$')
plt.semilogy(do.grid.t[1:], do.other.scalar.energyloss_T_cold[:,0], 'r', lw=3, label='DREAM')
plt.xlabel('Time $t$ (s)')
plt.ylabel(r'$q_\parallel$ and $P_{\rm halo}$')
plt.legend()
plt.tight_layout()

Aperp = qPar[1:] / (do.other.scalar.energyloss_T_cold[:,0])
#print(Aperp)

#print(do.other.scalar.energyloss_T_cold[:,0])


###
# Try to evaluate energyloss_T_cold separately
elT = do.grid.VpVol_f[-1] * do.other.scalar.Wcold_Tcold_Ar[:,0] * T[1:]
print(do.other.scalar.Wcold_Tcold_Ar[:,0])
#print(elT)

plt.figure()
plt.plot(do.grid.t[1:], elT, label='Estimated loss')
plt.plot(do.grid.t[1:], do.other.scalar.energyloss_T_cold[:,0], label='energyloss T cold')
plt.legend()

plt.figure()
plt.plot(do.grid.t[1:], Aperp)
plt.title(r'$A_\perp$')

plt.show()
