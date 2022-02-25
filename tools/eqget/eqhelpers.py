# Various helper routines for eqget

import numpy as np


def parametrize_equilibrium(psi_apRp, ptBPHI, ptx, pty, Rp, Zp, *args, **kwargs):
    """
    Parametrizes the given numeric equilibrium using the
    DREAM shaping parameters:

      kappa   -- Elongation
      delta   -- Triangularity
      Delta   -- Shafranov shift
      GOverR0 -- Toroidal magnetic field function
      psi     -- Poloidal flux
    """
    n = psi_apRp.size

    kappa, delta, Delta, GOverR0 = [], [], [], []
    for i in range(n):
        a = ptx[0,i]

        Zind = np.argmax(pty[:,i])
        R_upper = ptx[Zind,i]+Rp

        k = (np.amax(pty[:,i])-np.amin(pty[:,i])) / (2*a)
        #d = (ptx[0,i]+Rp-R_upper) / a
        D = (np.amax(ptx[:,i])+np.amin(ptx[:,i]))/2
        d = (Rp+D-R_upper) / a
        G = (ptx[0,i]+Rp)*ptBPHI[0,i] / Rp

        if np.isinf(k):
            k = 0

        kappa.append(k)
        delta.append(d)
        Delta.append(D)
        GOverR0.append(G)

    return {
        'R0': Rp,
        'r': np.array(ptx[0,:]),
        'kappa': np.array(kappa),
        'delta': np.array(delta),
        'Delta': np.array(Delta),
        'GOverR0': np.array(GOverR0),
        'psi': psi_apRp * ptx[0,-1]
    }


