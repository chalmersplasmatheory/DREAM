import numpy as np
import scipy.constants
from scipy import interpolate
import os, sys
import QAS as Stellarator

import Exceptions

sys.path.append('$DREAMPATH/py/')
sys.path.append('$DREAMPATH/py')
sys.path.append('$DREAMPATH/build/dreampyface/cxx/')
sys.path.append('$DREAMPATH/')
import dreampyface

from DREAM import DREAMSettings, runiface, DREAMOutput

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.OhmicCurrent as OhmicCurrent
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.TransportSettings as Transport
import DREAM.Settings.TimeStepper as TimeStepper
import DREAM.Settings.Equations.BootstrapCurrent as BootstrapCurrent


NR = 21
N_SAVE_STEPS = 200

PMAX = 2.5
NP = 80

# Globals
NT_MAX = 1000000
DT0_MIN = 1e-13

NT_ADD = 500000
DT_DIV = 10

DBB0 = 0#3e-3

def kineticGrids(ds, mode, pMax=PMAX, Np=NP, nr=NR):
    match mode:
        case 'fluid':
            ds.runawaygrid.setEnabled(False)
            ds.hottailgrid.setEnabled(False)
        case 'isotropic':
            ds.hottailgrid.setEnabled(True)
            ds.hottailgrid.setNxi(1)
            ds.hottailgrid.setPmax(pMax)
            ds.hottailgrid.setNp(Np)
        case 'superthermal':
            ds.hottailgrid.setEnabled(True)
            ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax=2 / nr)
            ds.hottailgrid.setPmax(pMax)
            ds.hottailgrid.setNp(Np)
        case _:
            raise Exception("ERROR: Not a valid mode")

    return ds

def terminateTQ(sim):
    tmin = 1e-5
    f_stop = 1 / 1000 * 4
    Tcold = sim.unknowns.getData('T_cold')['x']
    if sim.getCurrentTime() > tmin and (np.max(Tcold, axis=0) > Tcold[-1,:]).all():
        Tmean = np.mean(Tcold[-1,:])
        maxT = Stellarator.T_initial
        if Tmean <= f_stop * maxT:
            return True
    return False


def basicSetting(ds, nr=NR, pMax=PMAX, Np=NP, isotropic=False, superthermal=False, activated=None, bootstrap=False, out_init=None):
    """
    Define basic setting that are invariant for different simulation cases

    :param ds: DREAM settings object
    :return:
    """
    if superthermal:
        nr = 7

    # Temperature
    rT0, T0 = Stellarator.getInitialTemperatureProfile(nr)  # Initial profile

    # Initial densities
    if activated == 'DT':
        nD = 0.5 * Stellarator.getInitialDensityProfile(nr)[1]
        nT = 0.5 * Stellarator.getInitialDensityProfile(nr)[1]
    else:
        nD = 1. * Stellarator.getInitialDensityProfile(nr)[1]

    # Other parameters
    tau0 = 0#1 / Stellarator.tau_w  # inverse wall time

    # Set tokamak specifications
    Stellarator.setRadialGrid(ds, nr=nr)

    # Add initial ions to plasma
    if activated == 'DT':
        ds.eqsys.n_i.addIon('D', n=nD, Z=1, Z0=1, T=T0, r=rT0, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
        ds.eqsys.n_i.addIon('T', n=nT, Z=1, Z0=1, T=T0, r=rT0, tritium=True, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
    else:
        ds.eqsys.n_i.addIon('D', Z=1, Z0=1, iontype=Ions.IONS_DYNAMIC, T=T0, n=nD * np.ones((rT0.size,)), r=rT0, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
        ds.eqsys.n_re.setTritium(False)

    # Set collision settings
    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER  
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED  
    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT  
    ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL  

    # Use Sauter formula for conductivity
    ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)  

    ds = kineticGrids(ds, 'fluid')
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL  

    # Set solver settings
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
    ds.solver.setMaxIterations(maxiter=500)

    if isotropic or superthermal:
        ds.solver.tolerance.set(reltol=2e-4)
    else:
        ds.solver.tolerance.set(reltol=2e-4)
    ds.solver.tolerance.set('I_wall', abstol=1e-6)
    ds.solver.tolerance.set('V_loop_w', abstol=1e-10)
    ds.solver.tolerance.set('psi_wall', abstol=1e-6)
    ds.solver.tolerance.set(unknown='n_re', reltol=2e-6, abstol=1e5)
    ds.solver.tolerance.set(unknown='j_re', reltol=2e-6, abstol=1e-5)

    nfree, rn0 = ds.eqsys.n_i.getFreeElectronDensity()

    # Set current and current density profiles (instead of cond sim)
    rj, j = Stellarator.getInitialCurrentDensityProfile(nr)
    #Ip = Stellarator.Ip

    # Initial simulations to let E-field and plasma current etc stabilize at correct values
    # only used for isotropic and superthermal simulation, and if a non-uniform radial profile
    # is used for the injected material, this is used for normalization of the density.
    ds.eqsys.T_cold.setPrescribedData(T0, radius=rT0)

    ds.eqsys.E_field.setType(EField.TYPE_PRESCRIBED_OHMIC_CURRENT)
    ds.eqsys.j_ohm.setCurrentProfile(j, radius=rj)#, Ip0=Ip)

    ds.timestep.setTmax(1)
    ds.timestep.setNt(1)

    ds.solver.tolerance.set('n_hot', reltol=2e-6, abstol=1e5)
    ds.solver.tolerance.set('j_hot', reltol=2e-6, abstol=1)

    do_init = runiface(ds, out_init, quiet=True)

    if isotropic or superthermal:
        if superthermal:
            ds = kineticGrids(ds, 'superthermal', pMax=pMax, Np=Np, nr=nr)
        else:
            ds = kineticGrids(ds, 'isotropic', pMax=pMax, Np=Np)

        ds.eqsys.j_ohm.setCorrectedConductivity(False)  

        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=0.9999 * nfree, rT0=rT0, T0=T0)

        # All electrons between pThreshold and pMax are counted as hot electrons
        ds.eqsys.f_hot.setHotRegionThreshold(5)  
        # Boundary condition on f at p = pMax (assume f(p>pMax) = 0)
        ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
        # Enable flux limiters
        ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)
        ds.eqsys.f_hot.enableIonJacobian(False)  
    else:
        ds.eqsys.j_ohm.setInitialProfile(j, radius=rj)#, Ip0=Ip)
        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree, rT0=rT0, T0=T0)

    if bootstrap:
        # Include bootstrap current
        ds.eqsys.j_bs.setMode(BootstrapCurrent.BOOTSTRAP_MODE_ENABLED)
        ds.eqsys.j_bs.setInitMode(BootstrapCurrent.BOOTSTRAP_INIT_MODE_TOTAL)

    # Set fluid RE generation
    if superthermal:
        ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.01)
    else:
        ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    if isotropic or superthermal:
        ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)
    else:
        ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)
        ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_ANALYTIC_ALT_PC)

    if activated == 'DT':
        if isotropic or superthermal:
            ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_KINETIC, photonFlux=1.4e+18, C1=1.627, C2=0.850, C3=0.038)
            ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_KINETIC)
            ds.other.include(['fluid', 'scalar', 'hottail/S_compton', 'hottail/S_tritium'])

        else:
            ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_FLUID, photonFlux=1.4e+18, C1=1.627, C2=0.850, C3=0.038)
            ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_FLUID)
            ds.other.include(['fluid', 'scalar'])
    elif activated == 'DD':
        ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_KINETIC, photonFlux=3.3e+15, C1=1.525, C2=0.919, C3=0.094)
        ds.other.include(['fluid', 'scalar', 'hottail/S_compton'])
    else:
        ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_NEGLECT)
        ds.other.include(['fluid', 'scalar'])

    ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=tau0, R0=Stellarator.R0)

    # Set to self consistent temperature evolution
    ds.eqsys.T_cold.setType(Temperature.TYPE_SELFCONSISTENT)
    if isotropic or superthermal:
        ds.eqsys.T_cold.setInitialProfile(1)
    else:
        ds.eqsys.T_cold.setInitialProfile(T0, radius=rT0)

    if isotropic or superthermal:
        ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)

        ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)

        ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL

        ds.solver.tolerance.set('j_hot', reltol=1e-2, abstol=1e-5)
        ds.solver.tolerance.set('f_hot', reltol=1e-2, abstol=1e5)

        ignorelist = ['n_i', 'N_i', 'W_i', 'T_cold', 'W_cold', 'n_hot', 'n_cold']

        ds.fromOutput(out_init, ignore=ignorelist)

    return ds, do_init

# Will probably want this as injected over time
def injectMaterial(ds, nD, nNe, cD=0, cNe=0, do_init=None):
    if nD > 0:
        r_nD = np.linspace(0, Stellarator.a, NR)
        if cD == 0:
            profile = np.ones(NR)
            n_nD = nD * profile
        else:
            profile = (1 + np.tanh(cD * ((r_nD / Stellarator.a) - .5)))
            f_prof = interpolate.interp1d(r_nD, profile)
            profile_int = f_prof(do_init.grid.r)
            n_nD = nD * profile * do_init.grid.integrate(1) / do_init.grid.integrate(profile_int)
        ds.eqsys.n_i.addIon('D2', Z=1, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_nD, r=r_nD, T=1, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
    if nNe > 0:
        r_nNe = np.linspace(0, Stellarator.a, NR)
        if cNe == 0:
            profile = np.ones(NR)
            n_nNe = nNe * profile
        else:
            profile = (1 + np.tanh(cNe * ((r_nNe / Stellarator.a) - .5)))
            f_prof = interpolate.interp1d(r_nNe, profile)
            profile_int = f_prof(do_init.grid.r)
            n_nNe = nNe * profile * do_init.grid.integrate(1) / do_init.grid.integrate(profile_int)
        ds.eqsys.n_i.addIon('Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_nNe, r=r_nNe, T=1)
    return ds

# Other transport needed, if anything
def magneticPertubations(ds, dBB0, tmax, nt, nr=NR, kinetic=False):
    q = 1
    t = np.linspace(0, tmax, nt)
    r = np.linspace(0, Stellarator.a, nr)

    profile = np.ones(nr)
    dBB = dBB0 * profile

    Drr = np.pi * Stellarator.R0 * q * scipy.constants.c * dBB ** 2
    Drr = np.tile(Drr, (nt, 1))
    ds.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
    ds.eqsys.n_re.transport.prescribeDiffusion(Drr, t=t, r=r)

    if kinetic:
        ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)
        ds.eqsys.f_hot.transport.setMagneticPerturbation(dBB=np.tile(dBB, (nt, 1)), r=r, t=t)

    if dBB0 == 0:
        dBB = 4e-4 * profile
        ds.eqsys.n_re.transport.type = Transport.TRANSPORT_NONE
    R0dBB = dBB
    ds.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)
    ds.eqsys.T_cold.transport.setMagneticPerturbation(dBB=np.tile(R0dBB, (nt, 1)), r=r, t=t)
    return ds

def evalTQtime(do):
    f_stop = 1 / 1000
    t = do.grid.t[:]
    Tmean = np.mean(do.eqsys.T_cold[:, :], axis=1)
    imaxTm = np.argmax(Tmean)
    maxT = Stellarator.T_initial
    i = imaxTm
    for t_f, T_f in zip(t[imaxTm:], Tmean[imaxTm:]):
        if T_f <= f_stop * maxT:
            return t_f, i
        else:
            i += 1
    return None, -1

def runSimulation(nD_inj, nNe_inj, cD_inj=0, cNe_inj=0, isotropic=False, superthermal=False, activated=None, bootstrap=False, outputDir='.'):
    # Time parameters
    tmax_TQ_temp = Stellarator.tmax_TQ_max  # Maximum time for temporary thermal quench simulation
    tmax_CQ = Stellarator.tmax_CQ  # Total simulation time
    nt_TQ = 10000  # Temporary thermal quench simulation number of time steps
    nt_CQ = 20000  # Current quench number of time steps

    if isotropic or superthermal:
        dt0 = 1e-13
        kinetic = True
    else:
        dt0 = 1e-10
        kinetic = False
    dtmax = 1e-6

    dBB0_TQ = DBB0
    dBB0_CQ = 0.0

    # DREAM SETTINGS
    ds = DREAMSettings()

    # Basic DREAM settings
    out_init = f'{outputDir}/output_TQ_init.h5'
    ds, do_init = basicSetting(ds, isotropic=isotropic, superthermal=superthermal, activated=activated, bootstrap=bootstrap, out_init=out_init)

    ds = injectMaterial(ds, nD_inj, nNe_inj, cD_inj, cNe_inj, do_init=do_init)

    # Magnetic pertubations during thermal quench
    ds = magneticPertubations(ds, dBB0_TQ, tmax_TQ_temp, nt_TQ, kinetic=kinetic)

    ds.timestep.setType(TimeStepper.TYPE_IONIZATION)
    ds.timestep.setNt(None)
    ds.timestep.setTmax(None)

    # Test thermal quench
    ds.timestep.setIonization(dt0=dt0, dtmax=dtmax, tmax=tmax_TQ_temp)
    ds.timestep.setTerminationFunction(terminateTQ)
    out_TQ = f'{outputDir}/output_TQ.h5'
    ds.output.setFilename(out_TQ)
    s = dreampyface.setup_simulation(ds)
    #ds.save('settings_TQ.h5')
    try:
        do_TQ = s.run()#runiface(ds, out_TQ)
    except:
        ds.timestep.setIonization(dt0=dt0 / 10, dtmax=dtmax / 10, tmax=tmax_TQ_temp)
        s = dreampyface.setup_simulation(ds)

        try:
            do_TQ = s.run()
        except:
            raise Exceptions.DiscardSimulation("TQ simulation crashed")
    '''
    t_TQ, it_TQ = evalTQtime(do_TQ)
    
    if t_TQ == None:
        do_TQ.close()
        raise Exceptions.TransportException('Thermal quench was not successful.')
    else:
        ds.fromOutput(out_TQ, timeindex=it_TQ)
    '''
    t_TQ = do_TQ.grid.t[-1]
    if t_TQ >= tmax_TQ_temp:
        do_TQ.close()
        raise Exceptions.TransportException('Thermal quench was not successful.')

    ds.fromOutput(out_TQ)

    if t_TQ > 5e-4:
        nt_CQ = int(nt_CQ / 2)
    if t_TQ < 5e-5:
        nt_CQ = int(nt_CQ * 2)

    # Magnetic pertubations during current quench
    ds = magneticPertubations(ds, dBB0_CQ, tmax_CQ, nt_CQ, kinetic=kinetic)

    ds.solver.tolerance.set(unknown='n_re', reltol=2e-6, abstol=1e-15)

    if activated == 'DT':
        if kinetic:
            ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_KINETIC, photonFlux=1.4e+14, C1=1.627, C2=0.850, C3=0.038)
        else:
            ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_FLUID, photonFlux=1.4e+14, C1=1.627, C2=0.850, C3=0.038)
    elif activated == 'DD':
        ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_KINETIC, photonFlux=3.3e+11, C1=1.525, C2=0.919, C3=0.094)

    # Current quench
    ds.timestep.setType(TimeStepper.TYPE_CONSTANT)
    ds.timestep.setDt(None)
    ds.timestep.setNt(nt_CQ)
    ds.timestep.setTmax(tmax_CQ)
    ds.timestep.setNumberOfSaveSteps(int(nt_CQ / 10))
    ds.timestep.setTerminationFunction(None)
    out_CQ = f'{outputDir}/output_CQ.h5'
    ds.output.setFilename(out_CQ)
    s = dreampyface.setup_simulation(ds)
    try:
        do_CQ = s.run()#runiface(ds, out_CQ)
    except:
        ds.timestep.setNt(nt_CQ * 2)
        s = dreampyface.setup_simulation(ds)
        try:
            do_CQ = s.run()
        except:
            ds.timestep.setNt(nt_CQ * 4)
            s = dreampyface.setup_simulation(ds)
            try:
                do_CQ = s.run()
            except:
                raise Exceptions.DiscardSimulation("CQ simulation crashed")

    if kinetic:
        os.remove(out_init)

    return do_TQ, do_CQ
#runSimulation(0, 0, isotropic=False, activated='DT', outputDir='.')
runSimulation(5e21, 1e18, isotropic=False, activated='DT', bootstrap=True, outputDir='.')