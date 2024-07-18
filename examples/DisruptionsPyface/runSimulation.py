import argparse
import numpy as np
from pathlib import Path
import argparse
import scipy.constants
import sys
sys.path.append('../../py')
sys.path.append('../../build/dreampyface/cxx/')
sys.path.append('../..//Desktop/DREAM/')
import dreampyface
import DREAM

from DREAM import DREAMSettings, DREAMOutput, DREAMException, DREAMQuitException, runiface

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
import DREAM.Settings.Equations.DistributionFunction as DistFunc

import ITER as Tokamak

# Output file name pattern ('_XXX.h5' is appended to this)
PATTERN = 'output/'

# Radial resolution
NR = 21

# Time resolution
DT0 = [0, 1e-11, 1e-13, 1e-13, 1e-11] # TODO: Update
DTMAX = 1e-5
TMAX_TQ = Tokamak.tmax_TQ_max
TMAX_CQ = Tokamak.tmax_CQ_max
NT_CQ   = 20000

# Momentum resolution
NXI_CYLINDRICAL = 15
PMAX_KIN = 2.5  # maximum momentum (in mc) of hottail grid
PSEP_KIN = 0.3 # high grid-density region for p<psep (for mode SUPERTHERMAL, psep=0 effectively)
NP_KIN1 = 60 # number of momentum grid points below psep
NP_KIN2 = 80 # number of momentum grid points between psep and pmax (full Np in ISOTROPIC & SUPERTHERMAL mode)

# Magnetic perturbation strength
DBB0_TQ = 5e-3
DBB0_CQ = 0

# Available simulation modes
MODE_FLUID        = 1
MODE_SUPERTHERMAL = 2
MODE_KINETIC      = 3
MODE_ISOTROPIC    = 4

INCLUDE_FLUID_HOTTAIL = True    # Include fluid generation rates for hot-tail (in MODE_FLUID)
KINETIC_IONIZATION = True       # Kinetic ionization in non-fluid mode
SAUTER_CONDUCTIVITY = True      # Use the Sauter conductivity formula
SOLVER_RELTOL = 2e-6            # Default relative tolerance for Newton solver
LINEAR_SOLVER = Solver.LINEAR_SOLVER_MKL    # Linear solver to use

EXTENSION = ''


def setname(mode, phase, toroidal=True, io='settings', pattern="settings/"): # TODO: Keep toroidal? Keep different scenarios?
    return outname(mode=mode, phase=phase, toroidal=toroidal, io=io, pattern=pattern)


def outname(mode, phase, toroidal=True, io='output', pattern=PATTERN):  # TODO: Keep toroidal? Keep different scenarios?
    """
    Returns the appropriate output name to use.
    """
    base = '{}{}_{}_{}'.format(PATTERN, modename(mode), io, phase)

    if not toroidal:
        base += '_cyl'

    filename = '{}{}.h5'.format(base, EXTENSION)
    p = Path(filename).parent.resolve()

    if not p.exists():
        p.mkdir(parents=True)

    return filename


def modename(i):
    if i == MODE_FLUID: return 'fluid'
    elif i == MODE_KINETIC: return 'kinetic'
    elif i == MODE_SUPERTHERMAL: return 'superthermal'
    elif i == MODE_ISOTROPIC: return 'isotropic'
    else: return ''


def magneticPertubations(ds, dBB0, tmax, nt, mode=MODE_ISOTROPIC, nr=NR, REgrid=False):
    q = 1
    t = np.linspace(0, tmax, nt)
    r = np.linspace(0, Tokamak.a, nr)

    profile = np.ones(nr)
    dBB = dBB0 * profile

    Drr = np.pi * Tokamak.R0 * q * scipy.constants.c * dBB ** 2
    Drr = np.tile(Drr, (nt, 1))
    ds.eqsys.n_re.transport.setBoundaryCondition(Transport.BC_F_0)
    ds.eqsys.n_re.transport.prescribeDiffusion(Drr, t=t, r=r)

    if REgrid: # TODO
        ds.eqsys.f_re.transport.setBoundaryCondition(Transport.BC_F_0)
        ds.eqsys.f_re.transport.setMagneticPerturbation(dBB=np.tile(dBB, (nt, 1)), r=r, t=t)
    if mode != MODE_FLUID:
        ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)
        ds.eqsys.f_hot.transport.setMagneticPerturbation(dBB=np.tile(dBB, (nt, 1)), r=r, t=t)

    if dBB0 == 0:
        dBB = 4e-4 * profile
        ds.eqsys.n_re.transport.type = Transport.TRANSPORT_NONE
    R0dBB = dBB
    ds.eqsys.T_cold.transport.setBoundaryCondition(Transport.BC_F_0)
    ds.eqsys.T_cold.transport.setMagneticPerturbation(dBB=np.tile(R0dBB, (nt, 1)), r=r, t=t)
    return ds


def MMI(ds, name, Z, n, c=0., nr=NR):
    if n > 0:
        r_MMI = np.linspace(0, Tokamak.a, nr)
        if c == 0:
            profile = np.ones(nr)
            n_MMI = n * profile
        else:
            integrator, r_integrator = Tokamak.getIntegrator(nr=NR)
            profile = (1 + np.tanh(c * ((r_integrator / Tokamak.a) - .5)))
            n_MMI = n * profile * integrator.grid.integrate(1) / integrator.grid.integrate(profile)

        if Z == 1:
            ds.eqsys.n_i.addIon(name, Z=Z, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_MMI, r=r_MMI, T=1, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
            return ds
        ds.eqsys.n_i.addIon(name, Z=Z, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_MMI, r=r_MMI, T=1)
        return ds

def terminateTQ(sim):
    tmin = DTMAX
    f_stop = 1 / 1000
    Tcold = sim.unknowns.getData('T_cold')['x']
    if sim.getCurrentTime() > tmin and (np.max(Tcold, axis=0) > Tcold[-1,:]).all():
        Tmean = np.mean(Tcold[-1,:])
        maxT = Tokamak.T_initial
        if Tmean <= f_stop * maxT:
            return True
    return False


VpVol = None
dr = None
GR0 = None
Bmin = None
FSA_R02OverR2 = None

def terminateCQ(sim):
    # Tolerande parameters
    I_tol = Tokamak.Ire_tol

    # Get plasma, ohmic and RE currents
    Iohm = (VpVol * dr * GR0 / Bmin * FSA_R02OverR2 / (2 * np.pi) * sim.unknowns.getData('j_ohm')['x'][:, :]).sum(axis=1)
    Ire = (VpVol * dr * GR0 / Bmin * FSA_R02OverR2 / (2 * np.pi) * sim.unknowns.getData('j_re')['x'][:, :]).sum(axis=1)
    Ip = sim.unknowns.getData('I_p')['x'][:, :]

    tmax = sim.getCurrentTime()
    t = np.linspace(0, tmax, Iohm.shape[0])
    iIohmMax = np.argmax(Iohm)

    # Check if Ohmic current has decayed enough to evaluate CQ time and if there is a significant RE current
    if Iohm[-1] < 0.2 * Tokamak.Ip and Ire[-1] > I_tol:
        # Evaluate CQ time
        i_20Ip = np.argmin(np.abs(Iohm[iIohmMax:] - 0.2 * Tokamak.Ip)) + iIohmMax
        i_80Ip = np.argmin(np.abs(Iohm[iIohmMax:] - 0.8 * Tokamak.Ip)) + iIohmMax
        tauCQ = (t[i_20Ip] - t[i_80Ip]) / 0.6

        # Find out if the maximum or the 95 % RE current has occurred (representative RE current)
        Irerepr = (Ire[-40:] > 0.95 * Ip[-40:]).all() or (np.max(Ire) > Ire[-40:]).all()

        # End simulation if the representative RE current has occurred, and
        #   - if it has been running for longer than double of the evaluated CQ time, or
        #   - the ohmic current has reached below half of the tolerated RE current
        if (tmax > 2 * tauCQ or Iohm[-1] < I_tol/2) and Irerepr:
            return True
        return False
    # End simulation if the plasma current is below half of the tolerated RE current
    return Ip[-1] < I_tol/2

def runTQ(nD_MMI=0., cD_MMI=0., nNe_MMI=0., cNe_MMI=0., nAr_MMI=0., cAr_MMI=0., mode=MODE_ISOTROPIC, activated=True, toroidal=True): # TODO: Keep toroidal? Keep different scenarios?
    """
    TODO: Describe this
    """
    lsetname = lambda phase: setname(mode=mode, phase=phase, toroidal=toroidal)
    loutname = lambda phase: outname(mode=mode, phase=phase, toroidal=toroidal)

    nr = NR
    if mode == MODE_SUPERTHERMAL or mode == MODE_KINETIC:
        nr = (NR - 1)/2 + 1

    rT0, T0 = Tokamak.getInitialTemperatureProfile(Tokamak.T_core, Tokamak.T_edge, nr)  # Initial profile
    rj, j, Ip = Tokamak.getInitialCurrentDensityProfile(Tokamak.j_scale, Tokamak.j_exp, Tokamak.Ip, nr=nr)

    ds = DREAMSettings()

    # Set solver settings
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(LINEAR_SOLVER)
    ds.solver.setMaxIterations(maxiter=500)
    ds.solver.tolerance.set(reltol=2e-6)
    ds.solver.tolerance.set('I_wall', abstol=1e-6)
    ds.solver.tolerance.set('V_loop_w', abstol=1e-10)
    ds.solver.tolerance.set('psi_wall', abstol=1e-6)
    ds.solver.tolerance.set('n_re', abstol=1e5)
    ds.solver.tolerance.set('j_re', abstol=1e-5)
    if mode != MODE_FLUID:
        ds.solver.tolerance.set('n_hot', reltol=2e-6, abstol=1e5)
        ds.solver.tolerance.set('j_hot', reltol=2e-6, abstol=1)

    # Set collision settings
    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONAL # TODO: ok?

    # Set tokamak specifications
    Tokamak.setRadialGrid(ds, nr=nr, toroidal=toroidal)



    # Initial densities
    if activated:
        nD = 0.5 * Tokamak.ne0
        nT = 0.5 * Tokamak.ne0
        ds.eqsys.n_i.addIon('D', n=nD, Z=1, Z0=1, T=T0, r=rT0, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
        ds.eqsys.n_i.addIon('T', n=nT * np.ones((rT0.size,)), Z=1, Z0=1, T=T0, r=rT0, tritium=True, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
    else:
        nD = 1. * Tokamak.ne0
        ds.eqsys.n_i.addIon('D', Z=1, Z0=1, iontype=Ions.IONS_DYNAMIC, T=T0, n=nD * np.ones((rT0.size,)), r=rT0, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)

    nfree, rn0 = ds.eqsys.n_i.getFreeElectronDensity()

    if SAUTER_CONDUCTIVITY:
        ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)

    if mode == MODE_FLUID:
        ds.eqsys.j_ohm.setInitialProfile(j, radius=rj, Ip0=Ip)
        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree, rT0=rT0, T0=T0)
    else:
        # Initial simulations to let E-field and plasma current etc stabilize at correct values
        ds.runawaygrid.setEnabled(False)
        ds.hottailgrid.setEnabled(False)

        ds.eqsys.T_cold.setPrescribedData(T0, radius=rT0)
        ds.eqsys.E_field.setType(EField.TYPE_PRESCRIBED_OHMIC_CURRENT)
        ds.eqsys.j_ohm.setCurrentProfile(j, radius=rj, Ip0=Ip)

        ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

        ds.timestep.setTmax(1)
        ds.timestep.setNt(1)

        out_init = loutname('init')
        do_init = runiface(ds, out_init, quiet=True)
        ds.fromOutput(out_init, ignore=['n_i', 'N_i', 'W_i', 'T_cold', 'W_cold', 'n_hot', 'n_cold'])
        do_init.close()

        mod = 1
        if mode == MODE_SUPERTHERMAL or mode == MODE_ISOTROPIC:
            mod = 0.9999
            ds.eqsys.j_ohm.setCorrectedConductivity(False)

        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0= mod * nfree, rT0=rT0, T0=T0)

        # All electrons between pThreshold and pMax are counted as hot electrons
        ds.eqsys.f_hot.setHotRegionThreshold(5)
        # Boundary condition on f at p = pMax (assume f(p>pMax) = 0)
        ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
        # Enable flux limiters
        ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)
        ds.eqsys.f_hot.enableIonJacobian(False)

        if KINETIC_IONIZATION:
            ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)

        if mode == MODE_ISOTROPIC or mode == MODE_SUPERTHERMAL:
                ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL


    ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=1 / Tokamak.tau_w, R0=Tokamak.R0)

    # Set to self consistent temperature evolution
    ds.eqsys.T_cold.setType(Temperature.TYPE_SELFCONSISTENT)
    if mode == MODE_ISOTROPIC or mode == MODE_SUPERTHERMAL:
        ds.eqsys.T_cold.setInitialProfile(1)
    elif mode == MODE_FLUID:
        ds.eqsys.T_cold.setInitialProfile(T0, radius=rT0)


    match modename(mode):
        case 'fluid':
            ds.hottailgrid.setEnabled(False)
        case 'isotropic':
            ds.hottailgrid.setEnabled(True)
            ds.hottailgrid.setNxi(1)
            ds.hottailgrid.setPmax(PMAX_KIN)
            ds.hottailgrid.setNp(NP_KIN2)
        case 'superthermal':
            ds.hottailgrid.setEnabled(True)
            if not toroidal:
                ds.hottailgrid.setNxi(NXI_CYLINDRICAL)
            else:
                ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax=2 / nr)
            ds.hottailgrid.setPmax(PMAX_KIN)
            ds.hottailgrid.setNp(NP_KIN2)
        case 'kinetic':
            ds.hottailgrid.setEnabled(True)
            if not toroidal:
                ds.hottailgrid.setNxi(NXI_CYLINDRICAL)
            else:
                ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax=2 / nr)
            ds.hottailgrid.setPmax(PMAX_KIN)
            ds.hottailgrid.setNp(NP_KIN1 + NP_KIN2)
            ds.hottailgrid.setBiuniformGrid(psep=PSEP_KIN, npsep=NP_KIN1)

    include = ['fluid', 'scalar']

    if mode == MODE_FLUID or mode == MODE_ISOTROPIC:
        ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    else:
        ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.01)

    if mode == MODE_FLUID:
        ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)
        if INCLUDE_FLUID_HOTTAIL:
            ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_ANALYTIC_ALT_PC)
        else:
            ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_DISABLED)
    else:
        ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)
        ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_DISABLED)

    if activated:
        ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_KINETIC, photonFlux=1e18, C1=Runaways.C1_COMPTON_MS2017, C2=Runaways.C2_COMPTON_MS2017, C3=Runaways.C3_COMPTON_MS2017)
        ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_KINETIC)
        include.append('hottail/S_compton')
        include.append('hottail/S_tritium')
    else:
        ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_NEGLECT)
        ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_NEGLECT)

    MMI(ds, 'D_MMI', 1, nD_MMI, c=cD_MMI, nr=NR)
    MMI(ds, 'Ne', 10, nNe_MMI, c=cNe_MMI, nr=NR)
    MMI(ds, 'Ar', 18, nAr_MMI, c=cAr_MMI, nr=NR)

    magneticPertubations(ds, DBB0_TQ, TMAX_TQ, 1000, mode=mode, nr=nr)

    # Timestep settings
    ds.timestep.setType(TimeStepper.TYPE_IONIZATION)
    ds.timestep.setNt(None)
    ds.timestep.setTmax(None)
    ds.timestep.setIonization(dt0=DT0[mode], dtmax=DTMAX, tmax=TMAX_TQ)
    ds.timestep.setTerminationFunction(terminateTQ)

    ds.output.setFilename(loutname('TQ'))
    #ds.save(lsetname('TQ'))
    s = dreampyface.setup_simulation(ds)
    do = s.run()

    return ds, do

def runCQ(ds_TQ, do_TQ, mode=MODE_ISOTROPIC, activated=True, toroidal=True): # TODO: Keep toroidal? Keep different scenarios?
    lsetname = lambda phase: setname(mode=mode, phase=phase, toroidal=toroidal)
    loutname = lambda phase: outname(mode=mode, phase=phase, toroidal=toroidal)

    nr = NR
    if mode == MODE_SUPERTHERMAL or mode == MODE_KINETIC:
        nr = (NR - 1) / 2 + 1

    ds = DREAMSettings(ds_TQ)

    ds.fromOutput(loutname('TQ'))

    ds.solver.tolerance.set(unknown='n_re', reltol=2e-6, abstol=1e-15)

    if activated:
        ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_KINETIC, photonFlux=1e14, C1=Runaways.C1_COMPTON_MS2017, C2=Runaways.C2_COMPTON_MS2017, C3=Runaways.C3_COMPTON_MS2017)

    magneticPertubations(ds, DBB0_CQ, TMAX_CQ, 1000, mode=mode, nr=nr)

    global VpVol, dr, GR0, Bmin, FSA_R02OverR2
    VpVol = do_TQ.grid.VpVol[:]
    dr = do_TQ.grid.dr[:]
    GR0 = do_TQ.grid.GR0[:]
    Bmin = do_TQ.grid.Bmin[:]
    FSA_R02OverR2 = do_TQ.grid.FSA_R02OverR2[:]

    # Timestep settings
    ds.timestep.setType(TimeStepper.TYPE_CONSTANT)
    ds.timestep.setDt(None)
    ds.timestep.setNt(NT_CQ)
    ds.timestep.setTmax(TMAX_CQ)
    ds.timestep.setTerminationFunction(terminateCQ)
    ds.output.setFilename(loutname('CQ'))
    #ds.save(lsetname('CQ'))
    s = dreampyface.setup_simulation(ds)
    do = s.run()

    return ds, do

def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('-nD', help="MMI parameter: mean deuterium density", dest="nD", action='store', default=0., type=float)
    parser.add_argument('-cD', help="MMI parameter: deuterium density profile parameter", dest="cD", action='store', default=0., type=float)
    parser.add_argument('-nNe', help="MMI parameter: mean neon density", dest="nNe", action='store', default=0., type=float)
    parser.add_argument('-cNe', help="MMI parameter: deuterium neon profile parameter", dest="cNe", action='store', default=0., type=float)
    parser.add_argument('-nAr', help="MMI parameter: mean deuterium density", dest="nAr", action='store', default=0., type=float)
    parser.add_argument('-cAr', help="MMI parameter: deuterium density profile parameter", dest="cAr", action='store', default=0., type=float)

    parser.add_argument('-F', '--fluid', help="Use fluid mode", action='store_true', dest="fluid", default=False)
    parser.add_argument('-I', '--isotropic', help="Use isotropic mode", action='store_true', dest="isotropic", default=False)
    parser.add_argument('-S', '--superthermal', help="Use superthermal mode", action='store_true', dest="superthermal", default=False)
    parser.add_argument('-K', '--kinetic', help="Use kinetic mode", action='store_true', dest="kinetic", default=False)

    parser.add_argument('-A', '--activated', help="Activated plasma scenario", action='store_true', dest="activated", default=False)
    parser.add_argument('-C', '--cylindrical', help="Use cylindrical geometry", action='store_true', dest="cylindrical", default=False)
    settings = parser.parse_args()

    noMode = not np.array([settings.fluid, settings.isotropic, settings.superthermal, settings.kinetic]).any()

    if settings.nD == 0. and settings.nNe == 0. and settings.nAr == 0.:
        nD  = 1e22
        nNe = 1e19
    else:
        nD  = settings.nD
        nNe = settings.nNe
    nAr = settings.nAr

    if settings.cD >= -1 and settings.cD <= 1:
        cD  = settings.cD
    else:
        cD = 0.
    if settings.cNe >= -1 and settings.cNe <= 1:
        cNe = settings.cNe
    else:
        cNe = 0.
    if settings.cAr >= -1 and settings.cAr <= 1:
        cAr  = settings.cAr
    else:
        cAr = 0.

    if settings.fluid:
        print('Running in fluid mode')
        #ds_TQ, do_TQ = runTQ(nD_MMI=nD, cD_MMI=cD, nNe_MMI=nNe, cNe_MMI=cNe, nAr_MMI=nAr, cAr_MMI=cAr, mode=MODE_FLUID, activated=settings.activated, toroidal=(not settings.cylindrical))
        #ds_CQ, do_CQ = runCQ(ds_TQ, do_TQ, mode=MODE_FLUID, activated=settings.activated, toroidal=(not settings.cylindrical))

    if settings.isotropic or noMode:
        print('Running in isotropic mode')
        ds_TQ, do_TQ = runTQ(nD_MMI=nD, cD_MMI=cD, nNe_MMI=nNe, cNe_MMI=cNe, nAr_MMI=nAr, cAr_MMI=cAr, mode=MODE_ISOTROPIC, activated=settings.activated, toroidal=(not settings.cylindrical))
        ds_CQ, do_CQ = runCQ(ds_TQ, do_TQ, mode=MODE_ISOTROPIC, activated=settings.activated, toroidal=(not settings.cylindrical))

    if settings.superthermal:
        print('Running in superthermal mode')
        #ds_TQ, do_TQ = runTQ(nD_MMI=nD, cD_MMI=cD, nNe_MMI=nNe, cNe_MMI=cNe, nAr_MMI=nAr, cAr_MMI=cAr, mode=MODE_SUPERTHERMAL, activated=settings.activated, toroidal=(not settings.cylindrical))
        #ds_CQ, do_CQ = runCQ(ds_TQ, do_TQ, mode=MODE_SUPERTHERMAL, activated=settings.activated, toroidal=(not settings.cylindrical))

    if settings.kinetic:
        print('Running in kinetic mode')
        #ds_TQ, do_TQ = runTQ(nD_MMI=nD, cD_MMI=cD, nNe_MMI=nNe, cNe_MMI=cNe, nAr_MMI=nAr, cAr_MMI=cAr, mode=MODE_KINETIC, activated=settings.activated, toroidal=(not settings.cylindrical))
        #ds_CQ, do_CQ = runCQ(ds_TQ, do_TQ, mode=MODE_KINETIC, activated=settings.activated, toroidal=(not settings.cylindrical))



if __name__ == '__main__':
    sys.exit(main(sys.argv[:]))