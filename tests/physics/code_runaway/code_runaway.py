# CODE RUNAWAY
#
# Calculates the Dreicer runaway generation rate using the 
# full linearized collision operator in DREAM and compares 
# to the same calculation made with CODE.
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
import DREAM.Formulas

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.RunawayElectrons as Runaways

# Number of time steps to take
nTimeSteps = 4

def gensettings(T, Z=1, E=2, n=5e19, yMax=20):
    """
    Generate appropriate DREAM settings.

    T:    Electron temperature.
    E:    Effective charge of plasma.
    E:    Electric field (in units of critical electric field).
    n:    Electron density.
    yMax: Maximum momentum (normalized to thermal momentum) on
          computational grid.
    """
    betaTh = DREAM.Formulas.getNormalizedThermalSpeed(T)
    pMax = yMax * betaTh
    Ec = DREAM.Formulas.getEc(T, n)

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_THERMAL

    ds.eqsys.E_field.setPrescribedData(E)
    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   # Imaginary ion with charge Z
    ds.eqsys.n_cold.setPrescribedData(n)
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=n, rT0=0, T0=T)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

    ds.hottailgrid.setNxi(40)
    ds.hottailgrid.setNp(700)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setWallRadius(0.1)
    ds.radialgrid.setNr(1)

    tMax0 = pMax*Ec / E
    ds.timestep.setTmax(.9*tMax0)
    ds.timestep.setNt(nTimeSteps)

    ds.other.include('fluid/runawayRate', 'fluid/gammaDreicer')

    """ 
    # Alternativ parameters that run faster, but requires MUMPS
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_QUICK)
    ds.hottailgrid.setNxi(28)
    ds.hottailgrid.setNp(350)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MUMPS)
    """

    return ds


def getConnorHastieRate(T,E,ne=5e19,Z=1):
    """
    Calculates the Connor-Hastie runaway rate.
    """
    rrCH = np.zeros(E.shape)
    for i in range(0, E.size):
        rrCH[i] = DREAM.Formulas.getConnorHastieRunawayRate(T=T, n=ne, Zeff=Z, E=E[i])

    return rrCH


def loadCODE(filename):
    """
    Load conductivity values from CODE output.
    """
    T, E, rr = (None,)*3

    with h5py.File(filename, 'r') as f:
        T  = f['T'][:]
        E  = f['E'][:]
        rr = f['runawayRate'][:]

    return T, E, rr


def runTE(T, E):
    """
    Run DREAM for the specified values of temperature and ion charge.
    """
    ds = gensettings(T=T, E=E)
    ds.save('settings.h5')

    do = DREAM.runiface(ds, quiet=True)

    rrFull = do.other.fluid.runawayRate[:,0]
    rr     = rrFull[-1]
    # Connor-Hastie runaway rate
    rrCH   = do.other.fluid.gammaDreicer[-1,0]

    return rr, rrFull, rrCH


def run(args):
    """
    Run the test.
    """
    global nTimeSteps

    # Tolerance to require for agreement with CODE
    TOLERANCE = 5e-2    # 5%
    success = True
    workdir = pathlib.Path(__file__).parent.absolute()

    T, E, CODErr = loadCODE('{}/CODE-rates.mat'.format(workdir))

    nElong = 30

    nt     = nTimeSteps
    nE, nT = T.shape
    rr     = np.zeros((nE, nT))
    rrCH   = np.zeros((nE, nT))
    rrFull = np.zeros((nE, nT, nt))
    rrCH_E = np.zeros((nE, nT, nElong))
    rrCH   = np.zeros((nE, nT, nElong))
    for i in range(0, nE):
        for j in range(0, nT):
            print('Checking T = {} eV, E = {:.4f} V/m... '.format(T[i,j], E[i,j]), end="")
            try:
                rr[i,j], rrFull[i,j,:], rrCH[i,j] = runTE(T[i,j], E[i,j])
            except Exception as e:
                print(e)
                rr[i,j], rrFull[i,j,:] = 0, 0
                return False

            # Compare runaway rates
            Delta = 0
            if CODErr[i,j] < 1:
                # Some of the runaway rates are so small that they're unreliable,
                # so just compare to zero
                Delta = np.abs(rr[i,j])
            else:
                Delta = np.abs(rr[i,j] / CODErr[i,j] - 1.0)

            print("Delta = {:.3f}%".format(Delta*100))
            if Delta > TOLERANCE:
                dreamtests.print_error("DREAM runaway rate deviates from CODE at T = {} eV, E = {}".format(T[i,j], E[i,j]))
                success = False

            # Calculate Connor-Hastie rate
            rrCH_E[i,j,:] = np.linspace(E[0,j], E[-1,j], nElong)
            rrCH[i,j,:]   = getConnorHastieRate(T[0,j], rrCH_E[i,j])

    
    # Save
    if args['save']:
        with h5py.File('{}/DREAM-rates.h5'.format(workdir), 'w') as f:
            f['T'] = T
            f['E'] = E
            f['rr'] = rr

    if args['plot']:
        cmap = GeriMap.get()

        # Compare runaway rates
        plt.figure(figsize=(5,3.33))
        legs = []
        legh = []
        hN = None
        for i in range(0, nT):
            clr = cmap(i/nT)

            plt.loglog(rrCH_E[0,i,:], rrCH[0,i,:], color=clr, linewidth=2)

            h,  = plt.loglog(E[:,i], CODErr[:,i], 'o', color=clr, fillstyle='none', markersize=14, markeredgewidth=2)
            hN, = plt.loglog(E[:,i], rr[:,i], 'x', color=clr, markersize=10, markeredgewidth=3)

            s = r'$T = {:.0f}\,\mathrm{{eV}}$'.format(T[0,i])

            #legs.append(s)
            #legh.append(h)

        #legs.append('$\mathrm{DREAM}$')
        #legh.append(hN)

        plt.xlabel(r'Electric field $E_\parallel$ (V/m)')
        plt.ylabel(r'Runaway rate $\gamma\ \mathrm{(s)}^{-1}$')
        #plt.legend(legh, legs)

        plt.ylim([1e-28, 1e28])

        plt.tight_layout()
        plt.show()

    if success:
        dreamtests.print_ok("All runaway rates match those of CODE.")

    return success


