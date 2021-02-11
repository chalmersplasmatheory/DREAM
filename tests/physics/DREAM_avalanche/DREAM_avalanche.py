# DREAM AVALANCHE
#
# Calculates the avalanche growth rate using the superthermal
# collision operator in DREAM and compares with reference data
# generated using DREAM (with higher resolution)
# for a range of deuterium (neutral and ionised), argon and neon 
# plasma compositions as well as a wide range of electric fields.
#
# Reference data generated in commit: 
#   1df5375a25360784cc713a6f335945684e0e69d7 
# using resolution parameters:
#   nTimeSteps = 10
#   pOverPc    = 20
#   Nxi        = 50
#   Np         = 100 
#
# With flags --plot and/or --verbose, also compares the kinetic 
# simulations with the analytic (fluid) formulas implemented in DREAM.
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
    Generate DREAM settings object.

    T:          Electron temperature. (enters calculation only via lnLambda)
    EOverEcTot: Electric field (in units of critical electric field).
    nD0:        Density of ionised hydrogen.
    nD1:        Density of neutral hydrogen.
    nAr:        Density of neutral argon.
    nNe:        Density of neutral neon.
    """
    ######################
    # PHYSICAL CONSTANTS #
    ######################
    c    = scipy.constants.c            # speed of light
    ec   = scipy.constants.e            # elementary charge
    eps0 = scipy.constants.epsilon_0    # vacuum permittivity
    me   = scipy.constants.m_e          # electron mass

    #########################
    # RESOLUTION PARAMETERS #
    #########################
    pOverPc = 20  # pMax / pc, with pc an estimate of the critical momentum
    Nxi = 15      # number of xi grid points
    Np  = 60      # number of momentum grid points
    tMaxToP = 30  # time for collisionless acceleration to p/mc=tMaxToP

    ################################
    # SIMULATION PLASMA PARAMETERS #
    ################################
    nTot = nD0+nD1+nAr*18+nNe*10    # total (free plus bound) electron density
    nFree = nD1                     # free electron density

    lnLambda = 14.9-0.5*np.log(nFree/1e20) + np.log(T/1e3)
    EcTot = nTot*lnLambda*(ec**3) / (4*np.pi*(eps0**2)*me*(c**2))
    E = EcTot * EOverEcTot
    
    # Set pMax to a multiple of the critical momentum (in the nonscreened limit)
    # up to a maximum value of pMaxMax
    pcTot = 1/np.sqrt(E/EcTot-1)
    pMax = pOverPc * pcTot

    pMaxMax = 10
    if pMax>pMaxMax:
        pMax=pMaxMax
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
    
    # initialize f_hot to something small but smooth in order for the 
    # advection interpolation coefficients to converge but n_hot be 
    # negligible compared to n_re(t=0)
    ds.eqsys.f_hot.setInitialProfiles(n0=0.01, T0=1e5) 
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
    ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)

    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.01)
    ds.eqsys.n_re.setEceff(Eceff=Runaways.COLLQTY_ECEFF_MODE_SIMPLE)
    ds.eqsys.n_re.setInitialProfile(density=1) # arbitrary initial value for n_re to seed the avalanche
    ds.eqsys.f_hot.enableIonJacobian(False)

    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1e-6)
    ds.radialgrid.setMinorRadius(0.1)
    ds.radialgrid.setWallRadius(0.1)
    ds.radialgrid.setNr(1)

    tMax = tMaxToP*me*c / ((E-EcTot)*ec)
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(nTimeSteps)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.tolerance.set(reltol=1e-4)
    ds.solver.setVerbose(True)

    ds.other.include('fluid')

    return ds


def runNE(args,EOverEcTot=None, nD0=1e20, nD1=0, nAr=0, nNe=0):
    """
    Run DREAM for the specified values of electric field and ion density.
    
    Returns the avalanche growth rate (both numerical and analytic).
    """

    ds = gensettings(EOverEcTot=EOverEcTot, nD0=nD0, nD1=nD1, nAr=nAr, nNe=nNe)

    #ds.save('settings_DREAM_avalanche.h5')
    do = DREAM.runiface(ds, 'output.h5', quiet=True)

    GammaNumFull = do.other.fluid.runawayRate[:,0] / do.eqsys.n_re[1:,0]
    GammaNum     = GammaNumFull[-1]
    
    pMax = do.grid.hottail.p1_f[-1]
    pCrit = do.other.fluid.pCrit[0,0]
    pMaxOverPCrit = pMax/pCrit
    pMaxOverPCritCutOff = 3
    if args['verbose']:
        print('pMax/pCrit = {:.2f} (pMax = {:.2f}, pCrit = {:.2f}).'.format(pMaxOverPCrit, pMax, pCrit))
        var = abs(GammaNumFull[-1]/GammaNumFull[-2] - 1)
        if var > 1e-2:
            print('WARNING: growth rate may not be converged in time.')
            print('Variation in last two time steps: {:.2f}%'.format(100*var))
            if args['plot']:
                plotDiagnostics(do, GammaNumFull)
    if pMaxOverPCrit < pMaxOverPCritCutOff:
        print('WARNING: pMax/pCrit smaller than {:.3f}'.format(pMaxOverPCritCutOff))
        print('pMax/pCrit = {:.3f}.'.format(pMaxOverPCrit))


    '''
    Run two trivial fluid simulations in order to generate  
    growth rate data with the two different avalanche formulas
    '''
    ds.hottailgrid.setEnabled(False)
    ds.timestep.setNt(2)
    ds.solver.setType(Solver.LINEAR_IMPLICIT)

    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_FLUID)
    do = DREAM.runiface(ds, 'output.h5', quiet=True)
    GammaAn1Full = do.other.fluid.GammaAva[:,0]
    GammaAn1     = GammaAn1Full[-1]
    
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    do = DREAM.runiface(ds, 'output.h5', quiet=True)
    GammaAn2Full = do.other.fluid.GammaAva[:,0]
    GammaAn2     = GammaAn2Full[-1]

    return GammaNum, GammaNumFull, GammaAn1, GammaAn2

def run(args):
    """
    Run the test.
    """
    global nTimeSteps

    # Tolerance required for agreement with DREAM data.
    # The simulations are slightly underrresolved to speed up the test,
    # so we allow a 5% margin. 
    # Resolution bottlenecks appear to be Nxi, Np and pOverPc. 
    TOLERANCE = 5e-2
    success = True

    # Path to reference data file
    workdir = pathlib.Path(__file__).parent.absolute()    
    filename = '{}/DREAM-avalanche-data.h5'.format(workdir)

    # Define electric fields and densities to scan over
    EOverEcs = np.array([5,30,100])
    nDs = np.array([1e19,1e21])
    nZs = np.array([1e18,1e20])
    nE  = EOverEcs.size
    nnD = nDs.size
    nnZ = nZs.size

    nt = nTimeSteps

    # Unless run with --save, load reference data to compare with 
    if args['save']:
        GammaRef = np.zeros((nE, nnD, nnD, nnZ, nnZ))
    else:
        with h5py.File(filename, 'r') as f:
            GammaRef = f['GammaNum'][:]
            nDsRef = f['nD'][:]
            nZsRef = f['nZ'][:]
            EOverEcsRef = f['EOverEcTot'][:]

        # The plasma parameters must equal the saved reference data
        dataIsConsistent = np.array_equal(nDs,nDsRef) and np.array_equal(nZs,nZsRef) and np.array_equal(EOverEcs,EOverEcsRef)
        if not dataIsConsistent:
            dreamtests.print_error("Reference plasma parameters does not equal test parameters.")

    #                         E,  D0,  D1,  Ar,  Ne,  t
    GammaNum     = np.zeros((nE, nnD, nnD, nnZ, nnZ))
    GammaNumFull = np.zeros((nE, nnD, nnD, nnZ, nnZ, nt))
    GammaAn1     = np.zeros((nE, nnD, nnD, nnZ, nnZ))
    GammaAn2     = np.zeros((nE, nnD, nnD, nnZ, nnZ))
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
                            GammaNum[i,j,k,m,n], GammaNumFull[i,j,k,m,n:], GammaAn1[i,j,k,m,n], GammaAn2[i,j,k,m,n] = runNE(args,EOverEcTot=EOverEc, nD0=nD0, nD1=nD1, nAr=nAr, nNe=nNe)
                        except Exception as e:
                            print(e)
                            traceback.print_exc()
                            GammaNum[i,j,k,m,n], GammaNumFull[i,j,k,m,n:], GammaAn1[i,j,k,m,n], GammaAn2[i,j,k,m,n] = 0, 0, 0, 0
                            return False
                        if args['verbose']:
                            print('GammaNum = {:.3f}, GammaAva = {:.3f}, GammaAvaAlt = {:.3f}'.format(GammaNum[i,j,k,m,n],GammaAn1[i,j,k,m,n],GammaAn2[i,j,k,m,n]))
                        
                        # Compare growth rates
                        Delta    = np.abs( GammaNum[i,j,k,m,n] / GammaAn1[i,j,k,m,n] - 1.0)
                        DeltaAlt = np.abs( GammaNum[i,j,k,m,n] / GammaAn2[i,j,k,m,n] - 1.0)
                        
                        if args['save']: 
                            if args['verbose']:
                                print("Delta = {:.2f}%, DeltaAlt = {:.2f}%".format(Delta*100, DeltaAlt*100))
                        else:
                            DeltaRef = np.abs( GammaNum[i,j,k,m,n] / GammaRef[i,j,k,m,n] - 1.0)
                            if DeltaRef > TOLERANCE:
                                dreamtests.print_error("DREAM kinetic avalanche growth rate deviates from DREAM reference data")
                                success = False
                            if args['verbose']:
                                print("DeltaRef = {:.2f}%, Delta = {:.2f}%, DeltaAlt = {:.2f}%".format(DeltaRef*100, Delta*100, DeltaAlt*100))                        
                            else:
                                print("DeltaRef = {:.2f}%".format(DeltaRef*100))
                            
    # Plot or save results if requested
    if args['plot']:
        plotResults(GammaNum,GammaAn1,GammaAn2,EOverEcs, nDs, nZs, nE, nnD, nnZ)
    if args['save']:
        with h5py.File(filename, 'w') as f:
            f['GammaNum'] = GammaNum
            f['EOverEcTot'] = EOverEcs
            f['nD'] = nDs
            f['nZ'] = nZs
            dreamtests.print_ok("Test data saved to file as the future reference.")
    else:
        if success:
            dreamtests.print_ok("DREAM growth rate matches DREAM reference data.")

    return success

def plotResults(GammaNum,GammaAn1,GammaAn2,EOverEcs, nDs, nZs, nE, nnD, nnZ):
    """
    Visualizes the growth rates resulting from the test,
    comparing with analytic formulas
    """

    # Figure 1: plot growth rates vs E/Ec
    fig, axs = plt.subplots(2,3,num=1)

    Low=0
    High=-1
    plotSubplot(axs[0,0],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=High,iD1=Low,iAr=Low,iNe=Low, setLeg=True, setYLabel=True, fig=fig)
    axs[0,0].set_title(r'$n_\mathrm{{D}}^+ = {}$, others low'.format(nDs[High]))
    plotSubplot(axs[0,1],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=High,iD1=Low,iAr=High,iNe=Low)
    axs[0,1].set_title(r'$n_\mathrm{{D}}^+ = {}$, $n_\mathrm{{Ar}} = {}$'.format(nDs[High],nZs[High]))
    plotSubplot(axs[0,2],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=High,iD1=Low,iAr=Low,iNe=High)
    axs[0,2].set_title(r'$n_\mathrm{{D}}^+ = {}$, $n_\mathrm{{Ne}} = {}$'.format(nDs[High],nZs[High]))
    plotSubplot(axs[1,0],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=Low,iD1=High,iAr=Low,iNe=Low, setYLabel=True, setXLabel=True)
    axs[1,0].set_title(r'$n_\mathrm{{D}}^0 = {}$, others low'.format(nDs[High]))
    plotSubplot(axs[1,1],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=Low,iD1=High,iAr=High,iNe=Low, setXLabel=True)
    axs[1,1].set_title(r'$n_\mathrm{{D}}^0 = {}$, $n_\mathrm{{Ar}} = {}$'.format(nDs[High],nZs[High]))
    plotSubplot(axs[1,2],EOverEcs,GammaNum,GammaAn1,GammaAn2, iD0=Low,iD1=High,iAr=Low,iNe=High, setXLabel=True)#, setLeg=True, fig=fig)
    axs[1,2].set_title(r'$n_\mathrm{{D}}^0 = {}$, $n_\mathrm{{Ne}} = {}$'.format(nDs[High],nZs[High]))

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

    # Figure 2: scatter plot with Gamma_kinetic vs Gamma_fluid
    plt.figure(num=2)
    plotScatter(plt.gca(),GammaNum,GammaAn1,GammaAn2,nE, nnD, nnZ)

    plt.show()

def plotSubplot(ax,EOverEcs, GammaNum,GammaAn1,GammaAn2, iD0, iD1, iAr, iNe, setLeg=False, setXLabel=False, setYLabel=False, fig=None):
    """
    Plots a single subplot in the Gamma vs E/Ec results plot
    """

    l1,=ax.plot(EOverEcs,GammaNum[:,iD0,iD1,iAr,iNe], 'b' )
    l2,=ax.plot(EOverEcs,GammaAn1[:,iD0,iD1,iAr,iNe], 'r' )
    l3,=ax.plot(EOverEcs,GammaAn2[:,iD0,iD1,iAr,iNe], 'g' )

    if setXLabel:
        ax.set_xlabel(r'$E/E_{c,\mathrm{tot}}$')
    if setYLabel:
        ax.set_ylabel(r'$\Gamma$ [s$^{-1}$]')

    if setLeg and fig:
        ax.legend([l1,l2,l3],['DREAM kinetic','DREAM formula','NF 2019 formula'], loc="upper left")

def plotScatter(ax,GammaNum,GammaAn1,GammaAn2,nE, nnD, nnZ):
    """
    Plots a scatter plot of Gamma (numerical) vs Gamma (analytic) 
    """

    nLong = nE*nnD*nnD*nnZ*nnZ
    GammaNumLong = np.zeros(nLong)
    GammaAn1Long = np.zeros(nLong)
    GammaAn2Long = np.zeros(nLong)
    
    rms_An1 = 0
    rms_An2 = 0
    
    count=0
    for i in range(0, nE):
        for j in range(0, nnD):
            for k in range(0, nnD):
                for m in range(0, nnZ):
                    for n in range(0, nnZ):
                        Gn = GammaNum[i,j,k,m,n]
                        G1 = GammaAn1[i,j,k,m,n]
                        G2 = GammaAn2[i,j,k,m,n]
                        GammaNumLong[count] = Gn
                        GammaAn1Long[count] = G1
                        GammaAn2Long[count] = G2
                        rms_An1 = rms_An1 + (1-G1/Gn)**2
                        rms_An2 = rms_An2 + (1-G2/Gn)**2                        
                        count = count+1

    # Root-mean-square of the relative errors
    rms_An1 = np.sqrt(rms_An1/nLong)
    rms_An2 = np.sqrt(rms_An2/nLong)
    
    l1,=ax.plot(GammaNumLong[:],GammaAn1Long[:],'ro')
    l2,=ax.plot(GammaNumLong[:],GammaAn2Long[:],'go')
    
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    z1 = max( x1,y1 )
    z2 = min( x2,y2 )
    ax.plot([z1,z2],[z1,z2],'k--')

    t0 = ax.text(0.05,0.8, "DREAM RMS error: {:.2f}\% \n NF 2019 RMS error: {:.2f}\%".format(rms_An1*100,rms_An2*100),transform=ax.transAxes)
    t0.set_verticalalignment('top')
    t0.set_horizontalalignment('left')

    ax.set_xlabel(r'$\Gamma_\mathrm{kinetic}$')
    ax.set_ylabel(r'$\Gamma_\mathrm{formula}$')
    ax.legend((l1,l2),('DREAM','NF 2019'),loc='best')


def plotDiagnostics(do, GammaNumFull):
    """
    Plots diagnostic info about the solution 
    (mainly for debug)
    """
    print('pMax = {}, pCrit = {}'.format(do.grid.hottail.p1_f[-1], do.other.fluid.pCrit[0,0]))
    print('Ectot = {}, Eceff = {}'.format(do.other.fluid.Ectot[0,0],do.other.fluid.Eceff[0,0]))
    plt.figure(num=101)
    plt.plot(do.grid.t[1:],GammaNumFull)
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$\Gamma$ [s$^{-1}$]')

    plt.figure(num=102)
    plt.plot(do.grid.t,do.eqsys.n_re[:])
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$n_\mathrm{re}$ [m$^{-3}$]')

    plt.figure(num=105)
    mid_index = np.floor_divide(nTimeSteps,2)
    do.eqsys.f_hot.plot(t=[1,mid_index,-1],ax=plt.gca())

    plt.show()
