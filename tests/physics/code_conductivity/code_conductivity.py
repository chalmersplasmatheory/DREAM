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
import sys

import dreamtests

import DREAM
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMSettings import DREAMSettings
import DREAM.GeriMap as GeriMap

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.RunawayElectrons as Runaways


def gensettings(T, Z=1, EEc=1e-2, n=5e19, yMax=5):
    """
    Generate appropriate DREAM settings.

    T:    Electron temperature.
    Z:    Effective charge of plasma.
    EEc:  Electric field (in units of critical electric field).
    n:    Electron density.
    yMax: Maximum momentum (normalized to thermal momentum) on
          computational grid.
    """
    betaTh = DREAM.Formulas.getNormalizedThermalSpeed(T)
    pMax = yMax * betaTh
    Ec = DREAM.Formulas.getEc(T, n)

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_THERMAL

    ds.eqsys.E_field.setPrescribedData(EEc*Ec)
    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   # Imaginary ion with charge Z
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=n, rT0=0, T0=T)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

    ds.hottailgrid.setNxi(12)
    ds.hottailgrid.setNp(80)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setNr(1)

    # Simulate for 3.5 ms at T = 1 keV, and scale
    # appropriately depending on actual temperature
    tMax = 3.5e-3 * np.power(T / 1e3, 1.5)
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(2)
    
    return ds


def loadCODE(filename):
    """
    Load conductivity values from CODE output.
    """
    T, Z, sigma = (None,)*3

    with h5py.File(filename, 'r') as f:
        T     = f['T'][:]
        Z     = f['Z'][:]
        sigma = f['sigma'][:]

    return T, Z, sigma


def runTZ(T, Z):
    """
    Run DREAM for the specified values of temperature and ion charge.
    """
    ds = gensettings(T=T, Z=Z)
    E  = ds.eqsys.E_field[0,0]

    do = DREAM.runiface(ds, 'output.h5', quiet=True)
    j = do.eqsys.f_hot.currentDensity(t=-1)[0,0]
    sigma = j / E

    return sigma


def run(args):
    """
    Run the test.
    """
    # Tolerance to require for agreement with CODE
    TOLERANCE = 1e-2
    success = True
    workdir = pathlib.Path(__file__).parent.absolute()

    T, Z, CODEsigma = loadCODE('{}/CODE-conductivities.mat'.format(workdir))

    nZ, nT = T.shape
    sigma = np.zeros((nZ, nT))
    for i in range(0, nZ):
        for j in range(0, nT):
            print('Checking T = {} eV, Z = {:.1f}... '.format(T[i,j], Z[i,j]), end="")
            try:
                sigma[i,j] = runTZ(T[i,j], Z[i,j])
            except Exception as e:
                print('\x1B[1;31mFAIL\x1B[0m')
                print(e)
                sigma[i,j] = 0
                #continue
                return False

            # Compare conductivity right away
            Delta = np.abs(sigma[i,j] / CODEsigma[i,j] - 1.0)
            print("Delta = {:.3f}%".format(Delta*100))
            if Delta > TOLERANCE:
                dreamtests.print_error("DREAM conductivity deviates from CODE at T = {} eV, Z = {}".format(T[i,j], Z[i,j]))
                success = False
    
    # Save
    if args['save']:
        with h5py.File('{}/DREAM-conductivities.h5'.format(workdir), 'w') as f:
            f['T'] = T
            f['Z'] = Z
            f['sigma'] = sigma

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


