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

    # Poloidal flux is normalized by the aspect ratio, Rp/ap
    # (major radius over minor radius) in LUKE equilibria
    equil['psi_apRp'] = data['psi'] * ((data['R'][0,-1]-data['Rp']) / data['Rp'])
    equil['theta']    = data['theta']

    equil['ptx']      = data['R']-data['Rp']
    equil['pty']      = data['Z']-data['Zp']
    equil['ptBx']     = data['Br']
    equil['ptBy']     = data['Bz']
    equil['ptBPHI']   = data['Bphi']

    DREAMIO.SaveDictAsHDF5(filename, {'equil': equil})


