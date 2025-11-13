# Load equilibrium data from a file.

import numpy as np
import eqhelpers
import h5py

from DREAM.Settings.LUKEMagneticField import LUKEMagneticField
from EQDSK import EQDSK


def isAvailable():
    """
    Magnetic fields from file are always available, regardless of the
    computer system.
    """
    return True


def getLUKE(file, override_psilim=False, *args, **kwargs):
    """
    Returns magnetic equilibrium data from the named file.

    :param file: Name of file to load data from.
    """
    if h5py.is_hdf5(file):
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
    else:
        mf = EQDSK(file, override_psilim=override_psilim)
        return mf.get_LUKE(*args, **kwargs)


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


