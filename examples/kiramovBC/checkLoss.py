#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import e, atomic_mass

plt.rcParams.update({
    'font.family': 'FreeSans',
    'font.size': 15
})
import sys

sys.path.append('../../py/')
from DREAM import DREAMOutput


do = DREAMOutput('test.h5')

TIMESTEP = 10
kappa = 8
n = do.eqsys.n_cold[:,-1]
T = do.eqsys.T_cold[:,-1]
mi = 3.3435837724e-27
cs = np.sqrt(e * T / mi)
qPar = 2 * kappa * n[1:] * cs[1:] * e*T[1:]

LPar = 2 * np.pi * do.other.fluid.qR0[:,-1]

print(f'T = {e*T[TIMESTEP]}')
print(f'mi = {mi}')
print(f'Fr = {kappa * n[TIMESTEP] * cs[TIMESTEP] * e}')
#print(qPar)

S_wo_coeff = do.grid.VpVol_f[-1] / (do.grid.VpVol[-1] * do.grid.dr[-1])

plt.semilogy(do.grid.t[1:], do.other.scalar.Wcold_Tcold_Ar[:,0] * S_wo_coeff)
plt.semilogy(do.grid.t[1:], 4/3*kappa*n[1:]*cs[1:]*e / LPar, '--')


plt.figure(figsize=(6,4))

plt.semilogy(do.grid.t[1:], 2/3 * qPar / LPar, lw=3, label=r'$q_\parallel$')
plt.semilogy(do.grid.t[1:], do.other.scalar.energyloss_T_cold[:,0] / (do.grid.VpVol[-1] * do.grid.dr[-1]), 'r--', lw=3, label='DREAM')
plt.xlabel('Time $t$ (s)')
plt.ylabel(r'$q_\parallel$ and $P_{\rm halo}$')
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(do.grid.t[1:], do.other.scalar.energyloss_T_cold[:,0] / (do.grid.VpVol[-1] * do.grid.dr[-1]) / (2/3 * qPar / LPar))
plt.xlabel('Time $t$ (s)')
plt.ylabel(r'energylossTcold/qPar')

Aperp = qPar / (do.other.scalar.energyloss_T_cold[:,0])
#print(Aperp)

#print(do.other.scalar.energyloss_T_cold[:,0])


###
# Try to evaluate energyloss_T_cold separately
elT = do.grid.VpVol_f[-1] * do.other.scalar.Wcold_Tcold_Ar[:,0] * T[1:]
print(do.other.scalar.Wcold_Tcold_Ar[:,0])
#print(elT)

plt.figure()
plt.plot(do.grid.t[1:], elT, label='Estimated loss')
plt.plot(do.grid.t[1:], do.other.scalar.energyloss_T_cold[:,0], '--', label='energyloss T cold')
plt.legend()

plt.figure()
plt.plot(do.grid.t[1:], Aperp)
plt.title(r'$A_\perp$')

plt.show()
