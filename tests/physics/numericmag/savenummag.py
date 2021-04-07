# Save numerical magnetic field in given format


from DREAM import DREAMIO


def saveLUKE(filename, data):
    """
    Save the given numerical magnetic field in the LUKE format.

    :param data: Dictionary containing the raw magnetic field data.
    """
    equil = {'id': 'dream-numerical'}

    equil['Rp'] = data['Rp']
    equil['Zp'] = data['Zp']

    # TODO normalize
    equil['psi_apRp'] = data['psi']
    equil['theta']    = data['theta']

    equil['ptx']      = data['R']
    equil['pty']      = data['Z']
    equil['ptBx']     = data['Br']
    equil['ptBy']     = data['Bz']
    equil['ptBPHI']   = data['Bphi']

    DREAMIO.SaveDictAsHDF5(filename, {'equil': equil})


