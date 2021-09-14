# Load equilibrium data from a file.

import numpy as np
from DREAM.Settings.LUKEMagneticField import LUKEMagneticField


def isAvailable():
    """
    Magnetic fields from file are always available, regardless of the
    computer system.
    """
    return True


def getLUKE(file, *args, **kwargs):
    """
    Returns magnetic equilibrium data from the named file.

    :param file: Name of file to load data from.
    """
    mf = LUKEMagneticField(file)

    equil = {
        'id': mf.id,
        'Rp': mf.Rp,
        'Zp': mf.Zp,
        'psi_apRp': mf.psi_apRp,
        'theta': mf.theta,
        'ptx': mf.ptx,
        'pty': mf.pty,
        'ptBx': mf.ptBx,
        'ptBy': mf.ptBy,
        'ptBPHI': mf.ptBPHI
    }

    return equil

