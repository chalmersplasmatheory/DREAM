# CODE AVALANCHE
#
# Calculates the avalanche growth rate using the superthermal
# collision operator in DREAM and compares with analytic expressions
#  
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
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

# Number of time steps to take
nTimeSteps = 10
pOverPcTot = 3

def gensettings(T=10, EOverEcTot=None, nD0=1e20, nD1=0, nAr=0, nNe=0):
    """
    Generate appropriate DREAM settings.

    T:    Electron temperature.
    E:    Effective charge of plasma.
    E:    Electric field (in units of critical electric field).
    n:    Electron density.
    yMax: Maximum momentum (normalized to thermal momentum) on
          computational grid.
    """
    c    = scipy.constants.c
    ec   = scipy.constants.e
    eps0 = scipy.constants.epsilon_0
    me   = scipy.constants.m_e

    nTot = nD0+nD1+nAr*18+nNe*10
    nFree = nD1

    lnLambda = 14.9-0.5*np.log(nFree/1e20) + np.log(T/1e3)
    EcTot = nTot*lnLambda*(ec**3) / (4*np.pi*(eps0**2)*me*(c**2))
    E = EcTot * EOverEcTot
    
    # Set pMax to a multiple of the critical momentum (in the nonscreened limit)
    # up to a maximum value of pMaxMax
    pMaxMax = 2
    if E<(1/pMaxMax**2 + 1)*EcTot :
        pMax = pMaxMax
    else:
        pcTot = 1/np.sqrt(E/EcTot-1)
        pMax = pOverPcTot * pcTot

    ds = DREAMSettings()

    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

    ds.eqsys.E_field.setPrescribedData(E)

    ds.eqsys.n_i.addIon(name='D_ionized', Z=1, n=nD0, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED)   
    ds.eqsys.n_i.addIon(name='D_neutral', Z=1, n=nD1, iontype=IonSpecies.IONS_PRESCRIBED_NEUTRAL)   
    ds.eqsys.n_i.addIon(name='Ar', Z=18, n=nAr, iontype=IonSpecies.IONS_PRESCRIBED_NEUTRAL)   
    ds.eqsys.n_i.addIon(name='Ne', Z=10, n=nNe, iontype=IonSpecies.IONS_PRESCRIBED_NEUTRAL)   

    ds.eqsys.T_cold.setPrescribedData(T)
    ds.eqsys.f_hot.setInitialProfiles(n0=0, T0=T)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
    ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)

    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.01)
    ds.eqsys.n_re.setEceff(Eceff=Runaways.COLLQTY_ECEFF_MODE_SIMPLE)
    ds.eqsys.n_re.setInitialProfile(density=1)

    ds.hottailgrid.setNxi(20)
    ds.hottailgrid.setNp(30)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setNr(1)

    tMax = 20*me*c / (E*ec)
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(nTimeSteps)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setTolerance(1e-6)
    ds.solver.setVerbose(False)


    ds.other.include('fluid/runawayRate','fluid/GammaAva','fluid/GammaAvaAlt',)

    return ds


def runNE(EOverEcTot=None, nD0=1e20, nD1=0, nAr=0, nNe=0):
    """
    Run DREAM for the specified values of temperature and ion charge.
    """
    ds = gensettings(EOverEcTot=EOverEcTot, nD0=nD0, nD1=nD1, nAr=nAr, nNe=nNe)

    do = DREAM.runiface(ds, 'output.h5', quiet=True)

    GammaNumFull = do.other.fluid.runawayRate[:,0] / do.eqsys.n_re[1:,0]
    GammaNum     = GammaNumFull[-1]

    GammaAn1Full = do.other.fluid.GammaAva[:,0]
    GammaAn1     = GammaAn1Full[-1]
    
    GammaAn2Full = do.other.fluid.GammaAvaAlt[:,0]
    GammaAn2     = GammaAn2Full[-1]
    
    var = abs(GammaNumFull[-1]/GammaNumFull[-2] - 1)
    if var > 1e-2:
        print('Warning: growth rate not converged in time for')
        print('EOverEc = {}, nD0 = {} m-3, nD1 = {} m-3, nAr = {} m-3, nNe = {} m-3'.format(EOverEcTot, nD0, nD1, nAr, nNe))
        print('Variation in last two time steps of {}%'.format(100*var))
    return GammaNum, GammaNumFull, GammaAn1, GammaAn1Full, GammaAn2, GammaAn2Full


def run(args):
    """
    Run the test.
    """
    global nTimeSteps

    # Tolerance to require for agreement with analytic formula
    TOLERANCE = 20e-2    # 20%
    success = True

    # Define electric fields and densities to scan over
    nE  = 3
    '''
    nnD = 3
    nnZ = 3
    nDs = np.array([1e19,1e20,1e21])
    nZs = np.array([1e18,1e19,1e20])
    '''
    nnD = 2
    nnZ = 2
    nDs = np.array([1e19,1e21])
    nZs = np.array([1e18,1e20])
    EOverEcs = np.array([3,10,50])

    nt = nTimeSteps

    GammaNum     = np.zeros((nE, nnD, nnD, nnZ, nnZ))
    GammaNumFull = np.zeros((nE, nnD, nnD, nnZ, nnZ, nt))
    GammaAn1     = np.zeros((nE, nnD, nnD, nnZ, nnZ))
    GammaAn1Full = np.zeros((nE, nnD, nnD, nnZ, nnZ, nt))
    GammaAn2     = np.zeros((nE, nnD, nnD, nnZ, nnZ))
    GammaAn2Full = np.zeros((nE, nnD, nnD, nnZ, nnZ, nt))
    for i in range(0, nE):
        for j in range(0, nnD):
            for k in range(0, nnD):
                for m in range(0, nnZ):
                    for n in range(0, nnZ):
                        EOverEc=EOverEcs[i]
                        nD0=nDs[j]
                        nD1=nDs[k]
                        nAr=nZs[m]
                        nNe=nZs[n]
                        print('Checking E/Ectot = {}, nD0 = {} m-3, nD1 = {} m-3, nAr = {} m-3, nNe = {} m-3. '.format(EOverEc, nD0, nD1, nAr, nNe))
                        try:
                            GammaNum[i,j,k,m,n], GammaNumFull[i,j,k,m,n:], GammaAn1[i,j,k,m,n], GammaAn1Full[i,j,k,m,n:], GammaAn2[i,j,k,m,n], GammaAn2Full[i,j,k,m,n:] = runNE(EOverEcTot=EOverEc, nD0=nD0, nD1=nD1, nAr=nAr, nNe=nNe)
                        except Exception as e:
                            print(e)
                            GammaNum[i,j,k,m,n], GammaNumFull[i,j,k,m,n:], GammaAn1[i,j,k,m,n], GammaAn1Full[i,j,k,m,n:], GammaAn2[i,j,k,m,n], GammaAn2Full[i,j,k,m,n:] = 0, 0, 0, 0, 0, 0
                            return False
                        print('GammaNum = {}, GammaAva = {}, GammaAvaAlt = {}'.format(GammaNum[i,j,k,m,n],GammaAn1[i,j,k,m,n],GammaAn2[i,j,k,m,n]))
                        
                        # Compare growth rates
                        Delta    = np.abs( GammaNum[i,j,k,m,n] / GammaAn1[i,j,k,m,n] - 1.0)
                        DeltaAlt = np.abs( GammaNum[i,j,k,m,n] / GammaAn2[i,j,k,m,n] - 1.0)

                        print("Delta = {:f}%, DeltaAlt = {:f}%".format(Delta*100, DeltaAlt*100))
                        if Delta > TOLERANCE:
                            dreamtests.print_error("DREAM kinetic avalanche growth rate deviates from analytic formula")
                            success = False

    
    l1=plt.plot(GammaNum[:],GammaAn1[:])
    l2=plt.plot(GammaNum[:],GammaAn2[:])
    plt.xlabel('$\Gamma_\mathrm{num}$')
    plt.ylabel('$\Gamma_\mathrm{pred}$')
    plt.legend([l1,l2],('DREAM formula','Hesslow formula'))
    plt.show()
    # Save
    """
    with h5py.File('{}/DREAM-rates.h5'.format(workdir), 'w') as f:
        f['T'] = T
        f['E'] = E
        f['rr'] = rr
    """
    '''
    if args['plot']:
        cmap = GeriMap.get()

        # Compare runaway rates
        plt.figure(figsize=(9,6))
        legs = []
        legh = []
        hN = None
        for i in range(0, nT):
            clr = cmap(i/nT)

            h,  = plt.plot(E[:,i], CODErr[:,i], color=clr, linewidth=2)
            hN, = plt.plot(E[:,i], rr[:,i], 'x', color=clr, markersize=10, markeredgewidth=3)

            legs.append(r'$T = {:.0f}\,\mathrm{{eV}}$'.format(T[0,i]))
            legh.append(h)

        legs.append('$\mathrm{DREAM}$')
        legh.append(hN)

        plt.xlabel(r'$E$')
        plt.ylabel(r'$\sigma\ \mathrm{(S/m)}$')
        plt.legend(legh, legs)

        plt.show()
    '''
    if success:
        dreamtests.print_ok("All kinetic growth rates match the analytic formula.")

    return success


