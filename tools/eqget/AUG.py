
import h5py
import numpy as np
import sys
import eqhelpers

try:
    sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
    from sf2equ_20200525 import EQU
    import mapeq_20200507 as meq

    AVAILABLE = True
except:
    AVAILABLE = False


def isAvailable():
    """
    Returns ``True`` if this module can be used to fetch equilibrium data
    on this system.
    """
    global AVAILABLE
    return AVAILABLE


def getLUKE(shot, time, npsi=80, ntheta=80, filename=None):
    """
    Returns magnetic equilibrium data for the given time of the specified 
    AUG shot. If ``filename`` is provided, the data is also saved to the
    named LUKE equilibrium data file.

    The shape of the returned 2D arrays are (ntheta, npsi).

    :param shot: ASDEX Upgrade shot to fetch equilibrium data for.
    :param time: Time to fetch equilibrium data for.
    :param filename: Name of file to store data in.
    """
    equ = EQU(shot)

    # Radial grid (in normalized poloidal flux)
    rhop = np.linspace(0, 1, npsi+1)[1:]
    # Poloidal angle
    theta = np.linspace(0, 2*np.pi, ntheta)

    # Flux surface (R, Z) coordinates
    R, Z = meq.rhoTheta2rz(equ, rhop, theta, t_in=time, coord_in='rho_pol')
    R = R[0,:]
    Z = Z[0,:]

    # Poloidal flux psi
    psi = meq.rho2rho(equ, rhop, t_in=time, coord_in='rho_pol', coord_out='Psi')[0,:]

    # Calculate aspect ratio and normalize poloidal flux
    tidx = meq.get_nearest_index(equ.time, [time])[0][0]
    Rp   = equ.Rmag[tidx]
    Zp   = equ.Zmag[tidx]
    a    = R[0,-1]-Rp
    ieps = Rp / a

    psi_apRp = psi / ieps

    # Magnetic field components
    Br, Bz, Bphi = meq.rz2brzt(equ, r_in=R.flatten(), z_in=Z.flatten(), t_in=time)
    Br = Br[0,:].reshape(R.shape)
    Bz = Bz[0,:].reshape(R.shape)
    Bphi = Bphi[0,:].reshape(R.shape)

    equil = {
        'id': 'ASDEX Upgrade #{} t={:.4f}s'.format(shot, time),
        'Rp': np.array([Rp]), 'Zp': np.array([Zp]),
        'psi_apRp': psi_apRp,
        'theta': theta,
        'ptx': R-Rp, 'pty': Z-Zp,
        'ptBx': Br, 'ptBy': Bz, 'ptBPHI': Bphi
    }

    if filename:
        with h5py.File(filename, 'w') as f:
            f.create_group('equil')
            
            for key in equil.keys():
                f['equil/{}'.format(key)] = equil[key]

    return equil


def getShaping(shot, time, npsi=80, ntheta=80, filename=None, equil=None):
    if equil is None:
        equil = getLUKE(shot=shot, time=time, npsi=npsi, ntheta=ntheta)

    return eqhelpers.parametrize_equilibrium(**equil)


def getVolume(shot, time, filename=None):
    """
    Returns the plasma volume enclosed by a given flux surface.
    """
    tidx = meq.get_nearest_index(equ.time, [time])[0][0]

    data = {'psiN': equ.psiN[tidx,:], 'vol': equ.vol[tidx,:]}

    if filename:
        np.savez(filename, **data)

    return data
    
