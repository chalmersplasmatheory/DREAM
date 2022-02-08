# Load equilibrium data from a file.

import numpy as np
import eqhelpers

from DREAM.Settings.LUKEMagneticField import LUKEMagneticField
from GEQDSK import GEQDSK


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
    if file.endswith('.geqdsk'):
        mf = GEQDSK(file)
        return mf.get_LUKE(*args, **kwargs)
    else:
        # Assume LUKE format...
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

def getShaping(file, *args, equil=None, **kwargs):
    """
    Calculates DREAM shaping parameters corresponding to the given
    magnetic equilibrium.
    """
    if file.endswith('.geqdsk'):
        mf = GEQDSK(file)
        return mf.parametrize_equilibrium(*args, **kwargs)
    else:
        if equil is None:
            equil = getLUKE(file, *args, **kwargs)
        return eqhelpers.parametrize_equilibrium(**equil)


