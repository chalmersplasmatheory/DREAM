import numpy as np
import scipy.constants as const
from scipy.io import savemat
import sys, argparse, h5py
sys.path.append('../py')

from DREAM import DREAMOutput

def main():
    parser = argparse.ArgumentParser()

    # Arguments
    parser.add_argument('output', help="Name of DREAM output file (excluding .h5)", type=str)
    parser.add_argument('filename', help="Name of file to save data in (excluding file extension)", type=str)

    # Options
    parser.add_argument('-E', '--energy', help="Use energy pitch angle coordinates", action='store_const', const="E theta", dest="coordinates")
    parser.add_argument('-C', '--cylindrical', help="Use cylindrical coordinates", action='store_const',const="p_par p_perp", dest="coordinates")
    parser.set_defaults(coordinates="p xi")
    parser.add_argument('-M', '--matfile', help="Use Matlab version 5 file format to save data", action='store_true', dest="mat", default=False)
    parser.add_argument('-N', '--normalized', help="Use normalized distribution function", action='store_true', dest="normalized", default=False)
    parser.add_argument('-J', '--jacobian', help="Include jacobian in distribution function data to be saved", action='store_true', dest="jacobian", default=False)
    parser.add_argument('--time', help="Time indices for distribution function data.", nargs='*', type=int, default=[-1])

    settings = parser.parse_args()

    do = DREAMOutput(f'{settings.output}.h5')

    distFun = do.eqsys.f_re.volumeIntegratedDistributionFunction(settings.time, normalized=settings.normalized)

    t = do.grid.t[settings.time]
    p  = do.grid.runaway.P[:,:]
    xi = do.grid.runaway.XI[:,:]
    if settings.coordinates == 'p xi':

        if settings.jacobian:
            jac = (p**2).reshape(1,p.shape[0],p.shape[1])
            distFun *= jac

        dic = {"t": t, "p": p, "xi": xi, "f_re": distFun}

    elif settings.coordinates == 'E theta':

        mc2 = const.m_e * const.c**2
        e = const.e
        E = (np.sqrt((mc2*p)**2+mc2**2)-mc2) / e
        theta = np.arccos(xi)

        if settings.jacobian:
            dp_deE = e * (e*E + mc2) / (np.sqrt((e*E + mc2)**2 - mc2**2) * mc2)
            dxi_dtheta = np.sin(theta)
            jac = (p**2).reshape(1,p.shape[0],p.shape[1]) * dp_deE.reshape(1,p.shape[0],p.shape[1]) * dxi_dtheta.reshape(1,xi.shape[0],xi.shape[1])
            distFun *= jac

        dic = {"t": t, "E": E, "theta": theta, "f_re": distFun}

    elif settings.coordinates == 'p_par p_perp':

        p_par = do.grid.runaway.PPAR[:,:]
        p_perp = do.grid.runaway.PPERP[:,:]

        if settings.jacobian:
            jac = p**2 * (- p_perp / (p_par**2 + p_perp**2))
            distFun *= jac

        dic = {"t": t, "p_par": p_par, "p_perp": p_perp, "f_re": distFun}


    if settings.mat:
        savemat(f"{settings.filename}.mat", dic)
    else:
        with h5py.File(f"{settings.filename}.h5", 'w') as hf:
            for key, value in dic.items():
                if isinstance(value, (np.ndarray, list, tuple)):
                    hf.create_dataset("/" + key, data=np.array(value))
                else:
                    hf.create_dataset("/" + key, data=value)

    return 0

if __name__ == '__main__':
    sys.exit(main())