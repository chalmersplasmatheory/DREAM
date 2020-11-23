# CODE RUNAWAY
#
# Calculates the location in momentum p of the synchrotron bump as
# a function of magnetic field strength, B. This is done for a few
# different values of B and is then compared to the corresponding
# values calculated with CODE.
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
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

# Number of time steps to take
nTimeSteps = 4


def gensettings(B, T=5e3, Z=1, E=0.04, n=2e19, pMax=56):
    """
    Generate appropriate DREAM settings.

    B:    Magnetic field strength.
    T:    Electron temperature.
    E:    Effective charge of plasma.
    E:    Electric field (in units of critical electric field).
    n:    Electron density.
    yMax: Maximum momentum (normalized to thermal momentum) on
          computational grid.
    """
    global nTimeSteps

    c    = scipy.constants.c
    e    = scipy.constants.e
    eps0 = 8.85418782e-12
    me   = 9.10938e-31

    vth  = np.sqrt(2*e*T / me)

    lnLambda = 14.9-0.5*np.log(n/1e20) + np.log(T/1e3)
    Ec = n*lnLambda*(e**3) / (4*np.pi*(eps0**2)*me*(c**2))
    nu0 = n*e**4*lnLambda / (4*np.pi*eps0**2*me**2*(vth)**3)

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_THERMAL

    ds.eqsys.E_field.setPrescribedData(E)
    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   # Imaginary ion with charge Z
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
    ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)
    ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_QUICK)

    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

    Np = 105
    Nxi = 20
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)
    ds.hottailgrid.setBiuniformGrid(psep=1.5, npsep=20, thetasep=0.5, nthetasep_frac=0.5)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(B)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setNr(1)

    tMax_CODE = 600e3
    tMax = tMax_CODE / nu0
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(nTimeSteps)
    ds.output.setTiming(stdout=True)

    return ds


def loadCODE(filename):
    """
    Load conductivity values from CODE output.
    """
    B, bumpP = (None,)*2

    with h5py.File(filename, 'r') as f:
        B     = f['B'][:]
        bumpP = f['bumpP'][:]

    return B, bumpP


def findBump(do):
    # Locate bump
    tStart = 1
    t = do.grid.t[tStart:]
    p = np.zeros(t.shape)

    pCut = -19

    fPar = do.eqsys.f_hot[-1,0,-1,:pCut]

    df = np.diff(fPar)
    maxIdx = np.argmax(df)

    #plt.semilogy(do.grid.hottail.p[:pCut-1], np.abs(df))
    #plt.show()

    bumpLocIdx = np.argwhere(df[maxIdx:]<0)[0]
    return do.grid.hottail.p[maxIdx+bumpLocIdx]


def runB(B, pBump):
    """
    Run DREAM for the specified magnetic field strength.
    """
    ds = gensettings(B=B,pMax = min(56,2.5*pBump))
    do = DREAM.runiface(ds, quiet=True)

    # Locate synchrotron bump
    bumpP = findBump(do)

    return bumpP


def run(args):
    """
    Run the test.
    """
    global nTimeSteps

    # Tolerance to require for agreement with CODE
    TOLERANCE = 2e-2    # 2%
    success = True
    workdir = pathlib.Path(__file__).parent.absolute()

    B, CODEbump = loadCODE('{}/CODE-bumps.mat'.format(workdir))
    B = B[:,0]
    CODEbump = CODEbump[:,0]

    nt     = nTimeSteps
    nB     = B.size
    bumpP  = np.zeros((nB,))
    for i in range(1, nB):
    #for i in [2]:
        print('Checking B = {} T... '.format(B[i]), end="")
        try:
            bumpP[i] = runB(B[i],CODEbump[i])
        except Exception as e:
            print(e)
            bumpP[i] = 0
            return False

        # Compare runaway rates
        Delta = np.abs(bumpP[i] / CODEbump[i] - 1.0)

        print("Delta = {:.3f}%".format(Delta*100))
        if Delta > TOLERANCE:
            dreamtests.print_error("DREAM synchrotron bump location deviates from CODE at B = {} T".format(B[i]))
            dreamtests.print_error("DREAM bump at p = {} mc, CODE bump at p = {} mc".format(bumpP[i], CODEbump[i]))
            success = False

    
    # Save
    if args['save']:
        with h5py.File('{}/DREAM-bumps.h5'.format(workdir), 'w') as f:
            f['B'] = B
            f['bumpP'] = bumpP

    if args['plot']:
        plt.plot(B[1:], bumpP[1:], '*', markersize=10)
        plt.plot(B[1:], CODEbump[1:], 's', markersize=10)

        plt.legend(['DREAM', 'CODE'])
        plt.xlabel('Magnetic field strength (T)')
        plt.ylabel('Bump momentum ($mc$)')

        plt.tight_layout()
        plt.show()

    if success:
        dreamtests.print_ok("All runaway rates match those of CODE.")

    return success


