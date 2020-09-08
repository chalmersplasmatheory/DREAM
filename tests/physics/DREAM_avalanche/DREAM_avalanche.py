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
import traceback

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
nTimeSteps = 5

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
    ######################
    # PHYSICAL CONSTANTS #
    ######################
    c    = scipy.constants.c
    ec   = scipy.constants.e
    eps0 = scipy.constants.epsilon_0
    me   = scipy.constants.m_e

    #########################
    # RESOLUTION PARAMETERS #
    #########################
    pOverPc = 8   # pMax / pc, with pc an estimate of the critical momentum
    Nxi = 20      # number of xi grid points
    Np  = 40       # number of momentum grid points
    tMaxToP = 20  # time for collisionless acceleration to p/mc=tMaxInP

    ################################
    # SIMULATION PLASMA PARAMETERS #
    ################################
    nTot = nD0+nD1+nAr*18+nNe*10
    nFree = nD1

    lnLambda = 14.9-0.5*np.log(nFree/1e20) + np.log(T/1e3)
    EcTot = nTot*lnLambda*(ec**3) / (4*np.pi*(eps0**2)*me*(c**2))
    EcFree = nFree*lnLambda*(ec**3) / (4*np.pi*(eps0**2)*me*(c**2))
    E = EcTot * EOverEcTot
    
    # Set pMax to a multiple of the critical momentum (in the nonscreened limit)
    # up to a maximum value of pMaxMax
    pcTot = 1/np.sqrt(E/EcTot-1)
    pcFree = 1/np.sqrt(E/EcFree-1)
    pMax = pOverPc * np.sqrt(pcTot*pcFree) # set pc to the geometric mean of pcTot and pcFree
    if pMax>2:
        pMax=2
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
    ds.eqsys.f_hot.setInitialProfiles(n0=0, T0=T) # initialize f_hot = 0 
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
    ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)

    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.001)
    ds.eqsys.n_re.setEceff(Eceff=Runaways.COLLQTY_ECEFF_MODE_SIMPLE)
    ds.eqsys.n_re.setInitialProfile(density=1) # arbitrary initial value for n_re to seed the avalanche

    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setNr(1)

    tMax = tMaxToP*me*c / (E*ec)
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(nTimeSteps)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setTolerance(1e-4)
    ds.solver.setVerbose(True)
#    ds.solver.setType(Solver.LINEAR_IMPLICIT)

    ds.other.include('fluid')

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
    
    #plotDiagnostics(do, GammaNumFull)
    
    var = abs(GammaNumFull[-1]/GammaNumFull[-2] - 1)
    if var > 1e-2:
        print('WARNING: growth rate not converged in time for')
        print('EOverEc = {}, nD0 = {} m-3, nD1 = {} m-3, nAr = {} m-3, nNe = {} m-3'.format(EOverEcTot, nD0, nD1, nAr, nNe))
        print('Variation in last two time steps: {}%'.format(100*var))
        plotDiagnostics(do, GammaNumFull)
    return GammaNum, GammaNumFull, GammaAn1, GammaAn1Full, GammaAn2, GammaAn2Full

def plotDiagnostics(do, GammaNumFull):
    print('pMax = {}'.format(do.grid.hottail.p1[-1]))
    print('Ectot = {}, Eceff = {}'.format(do.other.fluid.Ectot[0],do.other.fluid.Eceff[0]))
    plt.figure(num=101)
    plt.plot(do.grid.t[1:],GammaNumFull)
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$\Gamma$ [s$^{-1}$]')

    plt.figure(num=102)
    plt.plot(do.grid.t,do.eqsys.n_re[:])
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$n_\mathrm{re}$ [m$^{-3}$]')

    plt.figure(num=103)
    do.eqsys.f_hot.plot(t=[1,np.round(.5*nTimeSteps),-1],ax=plt.gca())

    plt.show()


    

def run(args):
    """
    Run the test.
    """
    global nTimeSteps

    # Tolerance to require for agreement with analytic formula
    TOLERANCE = 20e-2    # 20%
    success = True

    # Define electric fields and densities to scan over
    nE  = 2
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
    EOverEcs = np.array([3,10])

    nt = nTimeSteps

    #                         E,  D0,  D1,  Ar,  Ne,  t
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
                            traceback.print_exc()
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

    plotResults(GammaNum,GammaAn1,GammaAn2,EOverEcs, nDs, nZs, nE, nnD, nnZ)

    if success:
        dreamtests.print_ok("All kinetic growth rates match the analytic formula.")

    return success

def plotResults(GammaNum,GammaAn1,GammaAn2,EOverEcs, nDs, nZs, nE, nnD, nnZ):
    
    # Figure 1: plot growth rates vs E/Ec
    fig, axs = plt.subplots(2,3,num=1)

    Low=0
    High=-1
    plotSubplot(axs[0,0],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=High,iD1=Low,iAr=Low,iNe=Low, setLeg=True, setYLabel=True, fig=fig)
    axs[0,0].set_title(r'$nD0 = {}$, others low'.format(nDs[High]))
    plotSubplot(axs[0,1],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=High,iD1=Low,iAr=High,iNe=Low)
    axs[0,1].set_title(r'$nD0 = {}$, $nAr = {}$'.format(nDs[High],nZs[High]))
    plotSubplot(axs[0,2],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=High,iD1=Low,iAr=Low,iNe=High)
    plotSubplot(axs[1,0],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=Low,iD1=High,iAr=Low,iNe=Low, setYLabel=True, setXLabel=True)
    plotSubplot(axs[1,1],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=Low,iD1=High,iAr=High,iNe=Low, setXLabel=True)
    plotSubplot(axs[1,2],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=Low,iD1=High,iAr=Low,iNe=High, setXLabel=True)#, setLeg=True, fig=fig)

    # Figure 2: scatter plot with Gamma_kinetic vs Gamma_fluid
#    ax = plt.figure(num=2)
#    plotScatter(ax,GammaNum,GammaAn1,GammaAn2,nE, nnD, nnZ)

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()

def plotSubplot(ax,EOverEcs, GammaNum,GammaAn1,GammaAn2, iD0, iD1, iAr, iNe, setLeg=False, setXLabel=False, setYLabel=False, fig=None):
    
    l1,=ax.plot(EOverEcs,GammaNum[:,iD0,iD1,iAr,iNe] )
    l2,=ax.plot(EOverEcs,GammaAn1[:,iD0,iD1,iAr,iNe] )
    l3,=ax.plot(EOverEcs,GammaAn2[:,iD0,iD1,iAr,iNe] )

    if setXLabel:
        ax.set_xlabel(r'$E/E_{c,\mathrm{tot}}$')
    if setYLabel:
        ax.set_ylabel(r'$\Gamma$ [s$^{-1}$]')

    if setLeg and fig:
        fig.legend([l1,l2,l3],['DREAM kinetic','DREAM formula','Hesslow formula'], loc="center right")

def plotScatter(ax,GammaNum,GammaAn1,GammaAn2,nE, nnD, nnZ):
    # do nothing
    nE=nE
    '''
    l1=ax.plot(GammaNumVec[:],GammaAn1Vec[:])
    l2=ax.plot(GammaNumVec[:],GammaAn2Vec[:])
    ax.xlabel(r'$\Gamma_\mathrm{num}$')
    ax.ylabel(r'$\Gamma_\mathrm{pred}$')
    ax.legend([l1,l2],('DREAM formula','Hesslow formula'))
    '''