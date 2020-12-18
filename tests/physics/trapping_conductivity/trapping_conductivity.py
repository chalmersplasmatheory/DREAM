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
import DREAM.Settings.Equations.OhmicCurrent as JOhm
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.RadialGrid as RadialGrid


NR = 3
R0 = 0.6

# in the chosen geometry and resolution, these are the trapped-passing boundaries
xi0Trapped = [0.39223227, 0.63245553, 0.76696499]

def gensettings(T, Z=300, EED=1e-6, n=5e19, yMax=5):
    """
    Generate appropriate DREAM settings.

    T:    Electron temperature.
    Z:    Effective charge of plasma.
    EEc:  Electric field (in units of critical electric field).
    n:    Electron density.
    yMax: Maximum momentum (normalized to thermal momentum) on
          computational grid.
    """
    global NR, R0

    betaTh = DREAM.Formulas.getNormalizedThermalSpeed(T)
    pMax = yMax * betaTh
    
    ED = DREAM.Formulas.getED(T,n)

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_THERMAL

    ds.eqsys.E_field.setPrescribedData(EED*ED)
    ds.eqsys.n_i.addIon(name='Ion', Z=Z, n=n/Z, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   # Imaginary ion with charge Z
    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=n, rT0=0, T0=T)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

    ds.eqsys.j_ohm.setCorrectedConductivity(JOhm.CORRECTED_CONDUCTIVITY_DISABLED)
    ds.eqsys.j_ohm.setConductivityMode(JOhm.CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS)

    # set non-uniform xi grid with cells stradding the trapped-passing boundaries
    ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(xi0Trapped, dxiMax=0.1, boundaryLayerWidth=1e-4)
    ds.hottailgrid.setNp(40)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    mu0 = scipy.constants.mu_0
    a  = 0.3
    rref = np.linspace(0, a, 20)
    Ipref = 2e5
    psiref = -mu0 * Ipref * (1-(rref/a)**2)

    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
    ds.radialgrid.setMajorRadius(R0)
    ds.radialgrid.setShaping(psi=psiref, rpsi=rref, kappa=1.5, delta=0.1, G=2.0)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a)
    ds.radialgrid.setNr(NR)

    # Simulate for 3.5 ms at T = 1 keV, and scale
    # appropriately depending on actual temperature
    tMax = 3.5e-3 * np.power(T / 1e3, 1.5)
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(3)
    
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


def runT(T):
    """
    Run DREAM for the specified values of temperature and ion charge.
    """
    global R0

    ds = gensettings(T=T, Z=300)
    ds.save('settings_trapping_conductivity.h5')
    do = DREAM.runiface(ds, quiet=True)
    jKinetic = do.eqsys.j_ohm[-1,:]

    ds.hottailgrid.setEnabled(False)
    do = DREAM.runiface(ds,quiet=True)
    jFluid = do.eqsys.j_ohm[-1,:]

    eps = do.grid.r[:] / R0

    return jKinetic, jFluid, eps


def run(args):
    """
    Run the test.
    """
    global NR

    # Tolerance to require for agreement with CODE
    TOLERANCE = 1e-2
    success = True
    workdir = pathlib.Path(__file__).parent.absolute()

    T = np.array([1e1,1e2,1e3,1e4])
    nT = T.size
    jKinetic = np.zeros((nT, NR))
    jFluid   = np.zeros((nT, NR))
    for i in range(0, nT):
        print('Checking T = {} eV... '.format(T[i]), end="")
        try:
            jKinetic[i,:], jFluid[i,:], eps = runT(T[i])
        except Exception as e:
            print('\x1B[1;31mFAIL\x1B[0m')
            print(e)
            #continue
            return False

        # Compare conductivity right away
        Delta = np.amax(np.abs(jKinetic[i,:] / jFluid[i,:] - 1.0))
        print("Delta = {:.3f}%".format(Delta*100))
        if Delta > TOLERANCE:
            dreamtests.print_error("DREAM conductivity deviates from CODE at T = {} eV".format(T[i]))
            success = False
    
    # Save
    if args['save']:
        with h5py.File('{}/DREAM-conductivities-trapping.h5'.format(workdir), 'w') as f:
            f['T'] = T
            f['jFluid']   = jFluid
            f['jKinetic'] = jKinetic

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

            h,  = plt.semilogy(eps, jKinetic[:,i], color=clr, linewidth=2)
            hN, = plt.semilogy(eps, jFluid[:,i], 'x', color=clr, markersize=10, markeredgewidth=3)

            if T[i] >= 1e3:
                s = r'$T = {:.0f}\,\mathrm{{keV}}$'.format(T[0,i]/1e3)
            else:
                s = r'$T = {:.0f}\,\mathrm{{eV}}$'.format(T[0,i])

            yfac = .6
            plt.text(0.5, jKinetic[0,i]*yfac, s, color=clr)

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


