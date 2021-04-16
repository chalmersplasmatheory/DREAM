#!/usr/bin/env python3
#
# This example generates a tokamak magnetic field using the analytica toroidal
# model, implemented separately in DREAM, and stores it as a numeric magnetic
# field. Two simulations are then run, utilizing both the analytic and
# corresponding numeric magnetic fields, and the results are compared.
#

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import scipy.constants
from scipy.interpolate import PchipInterpolator
import sys

import dreamtests
import numericmag.savenummag as savenummag

import DREAM
from DREAM import DREAMIO
from DREAM.DREAMOutput import DREAMOutput
from DREAM.DREAMSettings import DREAMSettings
import DREAM.GeriMap as GeriMap

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.RadialGrid as RadialGrid


# Path to parent directory of this script
ROOT = pathlib.Path(__file__).parent.absolute().resolve()


Rp = 2.0        # Major radius
a  = 0.55       # Minor radius
B0 = 6.2


def getShapeProfiles(nr=50):
    """
    Returns the shape profiles for magnetic field to use.
    """
    global a, B0, Rp

    r = np.linspace(0, a*1.05, nr)

    rG_R0, G_R0   = r, B0 * np.ones(r.shape)
    #rDelta, Delta = r, np.linspace(0, 0.05*a, nr)
    rkappa, kappa = r, np.linspace(1, 1.4, nr)
    rdelta, delta = r, np.linspace(0, 0.05, nr)
    rDelta, Delta = r, np.linspace(0, 0, nr)
    """
    rDelta, Delta = r, np.linspace(0, 0, nr)
    rkappa, kappa = r, np.linspace(1, 1, nr)
    rdelta, delta = r, np.linspace(0, 0, nr)
    """

    IpRef = 7e6
    mu0   = scipy.constants.mu_0
    rpsi  = r
    psi   = -mu0 * IpRef * (1-(rpsi/a)**2) * a

    return rG_R0, G_R0, rpsi, psi, rDelta, Delta, rkappa, kappa, rdelta, delta


def plotMagneticField(r, theta, R, Z, Br, Bz, Bphi, polar=False):
    """
    Plot the given magnetic field.
    """
    #fig, axs = plt.subplots(1,3, figsize=(24,8))
    fig, axs = plt.subplots(1,3, figsize=(12,4))

    yticks = np.pi * np.array([0, 0.5, 1, 1.5, 2])
    yticklabels = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

    if polar:
        x = r
        y = theta
    else:
        x = r*np.cos(theta)
        y = r*np.sin(theta)

    im = axs[0].contourf(x, y, Br)
    axs[0].set_title('Br')
    fig.colorbar(im, ax=axs[0])

    if polar:
        axs[0].set_yticks(yticks)
        axs[0].set_yticklabels(yticklabels)
    else:
        axs[0].set_aspect('equal', 'box')

    im = axs[1].contourf(x, y, Bz)
    axs[1].set_title('Bz')
    fig.colorbar(im, ax=axs[1])

    if polar:
        axs[1].set_yticks(yticks)
        axs[1].set_yticklabels(yticklabels)
    else:
        axs[1].set_aspect('equal', 'box')
    
    im = axs[2].contourf(x, y, Bphi)
    axs[2].set_title('Bphi')
    fig.colorbar(im, ax=axs[2])

    if polar:
        axs[2].set_yticks(yticks)
        axs[2].set_yticklabels(yticklabels)
    else:
        axs[2].set_aspect('equal', 'box')

    plt.show()


def constructMagneticField(Rp=2, Zp=0, a=0.5, nR=150, ntheta=151,
    rG_R0=None, G_R0=None, rpsi=None, psi=None,
    Delta=None, rDelta=None, kappa=None, rkappa=None,
    delta=None, rdelta=None, retdict=False):
    """
    Construct the numeric magnetic field.
    """
    raMax = 1.05
    r = np.linspace(0, a*1.05, nR+1)[1:]
    theta = np.linspace(0, 2*np.pi, ntheta+1)[:-1]

    mgR, mgT = np.meshgrid(r, theta)

    iG_R0, ipsi, iDelta, ikappa, idelta = (None,)*5

    if G_R0 is None: raise Exception('The toroidal magnetic field function must be specified.')
    else: iG_R0 = PchipInterpolator(rG_R0, G_R0, extrapolate=True)

    if psi is None: raise Exception('The poloidal flux must be specified.')
    else: ipsi = PchipInterpolator(rpsi, psi, extrapolate=True)

    if Delta is None: iDelta = PchipInterpolator([0, raMax], [0, 0], extrapolate=True)
    else: iDelta = PchipInterpolator(rDelta, Delta, extrapolate=True)

    if kappa is None: ikappa = PchipInterpolator([0, raMax], [1, 1], extrapolate=True)
    else: ikappa = PchipInterpolator(rkappa, kappa, extrapolate=True)

    if delta is None: idelta = PchipInterpolator([0, raMax], [0, 0], extrapolate=True)
    else: idelta = PchipInterpolator(rdelta, delta, extrapolate=True)

    R = lambda r, theta : Rp + iDelta(r) + r*np.cos(theta + idelta(r)*np.sin(theta))
    Z = lambda r, theta : Zp + r*ikappa(r)*np.sin(theta)

    psiPrime = ipsi.derivative()

    # Derivatives of shape parameters
    iDeltap = iDelta.derivative()
    ideltap = idelta.derivative()
    ikappap = ikappa.derivative()
    
    dRdr = lambda r, theta : iDeltap(r) + np.cos(theta + idelta(r)*np.sin(theta)) - r*ideltap(r)*np.sin(theta + idelta(r)*np.sin(theta))
    dZdr = lambda r, theta : ikappa(r) * (1 + r*ikappap(r)/ikappa(r)) * np.sin(theta)

    dRdt = lambda r, theta : -r*(1+idelta(r)*np.cos(theta)) * np.sin(theta + idelta(r)*np.sin(theta))
    dZdt = lambda r, theta : r*ikappa(r)*np.cos(theta)

    # Magnitude squared of the minor radius gradient
    def gradr2(r, theta):
        ct = np.cos(theta)
        st = np.sin(theta)
        sdt = np.sin(theta + idelta(r)*st)
        cdt = 1+idelta(r)*ct

        Jt1  = ikappa(r)*np.cos(idelta(r)*st)
        Jt2  = ikappa(r)*iDeltap(r)*ct 
        Jt3  = st*np.sin(theta+idelta(r)*st)
        Jt4  = r*ikappap(r) + ct*(idelta(r)*ikappa(r) + r*idelta(r)*ikappap(r) - r*ikappa(r)*ideltap(r))
        JOverRr = Jt1 + Jt2 + Jt3 * Jt4

        try:
            print('J/rR = {:.12f}'.format(JOverRr))
        except: pass
        return (ikappa(r)**2 * ct**2 + cdt**2 * sdt**2) / (JOverRr**2)

    Bphi = lambda r, theta : iG_R0(r) * Rp / R(r,theta)
    Br   = lambda r, theta : -psiPrime(r)*np.sqrt(gradr2(r,theta)) / (2*np.pi*R(r,theta)) * dZdr(r,theta) / np.sqrt(dRdr(r,theta)**2 + dZdr(r,theta)**2)
    Bz   = lambda r, theta : psiPrime(r)*np.sqrt(gradr2(r,theta)) / (2*np.pi*R(r,theta)) * dRdr(r,theta) / np.sqrt(dRdr(r,theta)**2 + dZdr(r,theta)**2)

    #B    = lambda r, theta : np.sqrt(Bphi(r,theta)**2 + Br(r,theta)**2 + Bz(r,theta)**2)

    if retdict:
        return {'Rp': Rp, 'Zp': Zp, 'psi': ipsi(r), 'theta': theta, 'R': R(mgR,mgT), 'Z': Z(mgR,mgT), 'Br': Br(mgR,mgT), 'Bz': Bz(mgR,mgT), 'Bphi': Bphi(mgR, mgT)}
    else:
        return Rp, Zp, ipsi(r), theta, R(mgR,mgT), Z(mgR,mgT), Br(mgR,mgT), Bz(mgR,mgT), Bphi(mgR,mgT)


def generateSettings(analyticB=False):
    """
    Generates a DREAMSettings object for the test

    :param bool analyticB: If ``True``, uses an analytic magnetic field. Otherwise a numeric magnetic field is used.
    """
    global ROOT

    T = 3e3     # eV
    E = 2       # V/m
    n = 5e19    # m^-3
    yMax = 20   # thermal momentum
    Z = 1       # plasma charge

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
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_QUICK)
    
    ds.hottailgrid.setNxi(20)
    ds.hottailgrid.setNp(100)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    # Get magnetic field shaping parameters
    rG_R0, G_R0, rpsi, psi, rDelta, Delta, rkappa, kappa, rdelta, delta = getShapeProfiles()
    if analyticB:
        ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
        ds.radialgrid.setShaping(psi=psi/Rp, rpsi=rpsi, GOverR0=G_R0, rG=rG_R0, kappa=kappa, rkappa=rkappa, delta=delta, rdelta=rdelta, Delta=Delta, rDelta=rDelta)
    else:
        FAC = 5
        numdata = constructMagneticField(Rp=Rp, Zp=0, a=a, nR=FAC*50, ntheta=FAC*50,
            rG_R0=rG_R0, G_R0=G_R0, rpsi=rpsi, psi=psi,
            Delta=Delta, rDelta=rDelta, kappa=kappa, rkappa=rkappa,
            delta=delta, rdelta=rdelta, retdict=True)

        # Save to HDF5 file
        numname = '{}/magfield.h5'.format(ROOT)
        savenummag.saveLUKE(numname, numdata)

        ds.radialgrid.setNumerical(numname, format=RadialGrid.FILE_FORMAT_LUKE)

    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a)
    ds.radialgrid.setMajorRadius(Rp)
    ds.radialgrid.setNr(10)

    tMax0 = pMax*Ec / E
    ds.timestep.setTmax(.9*tMax0)
    ds.timestep.setNt(5)

    ds.other.include('fluid/runawayRate', 'fluid/gammaDreicer')

    return ds


def runSimulation(analyticB):
    """
    Run a simulation and return the resulting value of the
    conductivity.
    """
    global ROOT

    t = 'analytical' if analyticB else 'numerical'

    ds = generateSettings(analyticB)
    ds.save('{}/settings_{}.h5'.format(ROOT, t))

    do = DREAM.runiface(ds, quiet=True)
    j  = do.eqsys.f_hot.currentDensity(t=-1)[0,:]
    sigma = j / do.eqsys.E_field[-1,:]

    return sigma


def run(args):
    """
    Run the test.
    """
    global ROOT

    TOLERANCE = 1e-2
    success = True

    ds = generateSettings(False)
    ds.save('{}/settings_numerical.h5'.format(ROOT))

    
    print('Comparing conductivity in analytical and numerical cases... ')

    print(' - Numerical... ')
    sigmaN = runSimulation(False)
    print(' - Analytical... ')
    sigmaA = runSimulation(True)

    Delta = np.max(np.abs(sigmaN / sigmaA - 1))
    print(" ==> Delta = {:.3f}%".format(Delta*100))
    if Delta > TOLERANCE:
        print(' \x1B[1;31mFAIL\x1b[0m')
        success = False
    

    return success

