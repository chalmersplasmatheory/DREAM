# CODE CONDUCTIVITY
#
# Calculates the conductivity using the full linearized collision operator in
# DREAM and compares to the same calculation made with CODE.
#
############################################################################

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import scipy.constants
import scipy.special
import sys

#import dreamtests

import DREAM
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMSettings import DREAMSettings
import DREAM.GeriMap as GeriMap
from DREAM.Formulas.PlasmaParameters import evaluateBraamsConductivity

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.RunawayElectrons as Runaways


def gensettings(T, Z=1, n=5e19):
    """
    Generate appropriate DREAM settings.

    T:    Electron temperature.
    Z:    Effective charge of plasma.
    n:    Electron density.
    """
    a = 0.5
    nr = 100

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_THERMAL

    # Get initial electric field profile
    r0, E0, lmbd = getInitialElectricField()
    ds.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setInitialProfile(efield=E0, radius=r0*a)
    ds.eqsys.E_field.setBoundaryCondition(bctype=Efield.BC_TYPE_SELFCONSISTENT, inverse_wall_time=0, R0=1)

    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   # Imaginary ion with charge Z
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.n_re.setDreicer(False)
    ds.eqsys.n_re.setAvalanche(False)

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a)
    ds.radialgrid.setNr(nr)

    # Expected e-folding time scale is mu0*sigma
    tMax = scipy.constants.mu_0 * evaluateBraamsConductivity(n, T, Z) / (lmbd[0]/a)**2 * 100
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(100)

    ds.other.include('fluid/conductivity')
    
    return ds, lmbd


def getInitialElectricField(nr=100, A=0.1):
    r"""
    Constructs the initial electric field profile

      E(r) = sum_k A_k*J_0(\lambda_k*r)

    from the given coefficients A_k. Here, J_0(x) is the zeroth
    Bessel function and \lambda_k is its k'th zero.

    :param float a: Plasma minor radius.
    :param float A: List of coefficients

    :return: Radial grid, initial electric field profile, and
    list of Bessel zeros (one per coefficient A_k specified).
    """
    if np.isscalar(A):
        A = np.array([A])

    lmbd = scipy.special.jn_zeros(0, A.size)
    
    r  = np.linspace(0, 1, nr)
    E0 = np.zeros((nr,))

    for k in range(0,A.size):
        E0 += A[k] * scipy.special.j0(lmbd[k]*r)

    return r, E0, lmbd


def runTZ(T, Z, full_output=False):
    """
    Run DREAM for the specified values of temperature and ion charge.
    """
    ds, lmbd = gensettings(T=T, Z=Z)
    do = DREAM.runiface(ds, 'output.h5', quiet=True)

    E0 = do.eqsys.E_field[0,:]
    exp = do.eqsys.E_field[1:,:] / (E0*scipy.special.j0(lmbd[0]*do.grid.r[:]/ds.radialgrid.a))

    if len(lmbd) != 1:
        raise Exception("Too many Bessel modes selected for the simulation. This test can only handle one mode.")

    t     = do.grid.t[1:].reshape((exp.shape[0], 1))
    #sigma = -scipy.constants.mu_0 * np.log(exp) / (lmbd[0]**2 * t)
    sigma = -(lmbd[0]/ds.radialgrid.a)**2 * t / (scipy.constants.mu_0 * np.log(exp))

    s      = do.other.fluid.conductivity[:]
    dSigma = np.abs((s - sigma) / s)

    plt.plot(t, sigma[:,-1])
    plt.plot(t, s[:,-1], 'r--')
    plt.show()

    if full_output:
        return do.grid.r, dSigma, sigma, s
    else:
        return dSigma


def run(args):
    """
    Run the test.
    """
    # Tolerance to require for agreement with CODE
    TOLERANCE = 1e-2
    success = True
    workdir = pathlib.Path(__file__).parent.absolute()

    #T, Z = np.meshgrid([100, 1000, 2000, 5000, 10000], [1, 2, 5, 10, 20, 50])
    T, Z = np.meshgrid([1000], [1])

    nZ, nT = T.shape
    exp = np.zeros((nZ, nT))
    for i in range(0, nZ):
        for j in range(0, nT):
            print('Checking T = {} eV, Z = {:.1f}... '.format(T[i,j], Z[i,j]), end="")
            try:
                exp[i,j] = runTZ(T[i,j], Z[i,j])
            except Exception as e:
                print('\x1B[1;31mFAIL\x1B[0m')
                print(e)
                exp[i,j] = 0
                #continue
                return False

            # Compare conductivity right away
            Delta = np.abs(exp[i,j] / CODEsigma[i,j] - 1.0)
            print("Delta = {:.3f}%".format(Delta*100))
            if Delta > TOLERANCE:
                dreamtests.print_error("DREAM conductivity deviates from CODE at T = {} eV, Z = {}".format(T[i,j], Z[i,j]))
                success = False
    
    # Save
    if args['save']:
        with h5py.File('{}/DREAM-conductivities.h5'.format(workdir), 'w') as f:
            f['T'] = T
            f['Z'] = Z
            f['exp'] = exp

    if args['plot']:
        cmap = GeriMap.get()

        # Compare conductivities
        #plt.figure(figsize=(9,6))
        plt.figure(figsize=(5,3))
        legs = []
        legh = []
        hN = None
        for i in range(0, nT):
            clr = cmap(i/nT)

            h,  = plt.semilogy(Z[:,i], CODEsigma[:,i], color=clr, linewidth=2)
            hN, = plt.semilogy(Z[:,i], sigma[:,i], 'x', color=clr, markersize=10, markeredgewidth=3)

            if T[0,i] >= 1e3:
                s = r'$T = {:.0f}\,\mathrm{{keV}}$'.format(T[0,i]/1e3)
            else:
                s = r'$T = {:.0f}\,\mathrm{{eV}}$'.format(T[0,i])

            yfac = 1.2
            plt.text(Z[-2,i]+4, sigma[-2,i]*yfac, s, color=clr)

            legs.append(s)
            legh.append(h)

        legs.append('$\mathrm{DREAM}$')
        legh.append(hN)

        plt.xlabel(r'$Z$')
        plt.ylabel(r'$\sigma\ \mathrm{(S/m)}$')
        #plt.legend(legh, legs)

        plt.tight_layout()
        plt.show()

    if success:
        dreamtests.print_ok("All conductivities match those of CODE.")

    return success


if __name__ == '__main__':
    r, dsigma, sigma_calc, sigma = runTZ(100, 2, full_output=True)

    #print('dSigma = {}'.format(dsigma[0]))
    #print(f'sigma  = {sigma_calc}')
    #print(f'braams = {sigma}')
    #plt.plot(r, sigma_calc[-1,:])
    #plt.plot(r, sigma[-1,:])
    plt.plot(r, sigma[-1,:] / sigma_calc[-1,:])
    print(sigma[-1,0] / sigma_calc[-1,0])
    #plt.legend([r'Calculated $\sigma$', r'Braams-Karney $\sigma$'])
    plt.show()

