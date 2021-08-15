#!/usr/bin/env python3

import argparse
import numpy as np
from pathlib import Path
import sys

try:
    import DREAM
except ModuleNotFoundError:
    sys.path.append('../../py')
    import DREAM

from DREAM import DREAMSettings

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.OhmicCurrent as OhmicCurrent
import DREAM.Settings.Equations.RunawayElectrons as RunawayElectrons
import DREAM.Settings.Solver as Solver

import ASDEXU as Tokamak
#import ITER as Tokamak

# Output file name pattern ('_XXX.h5' is appended to this)
PATTERN = 'output/out'

# Radial resolution
NR = 15

# Kinetic grid resolution
NP_KIN1 = 60 # number of momentum grid points below psep
NP_KIN2 = 80 # number of momentum grid points between psep and pmax (full Np in ISOTROPIC & SUPERTHERMAL mode)

# nXi used in cylindrical mode (in toroidal mode the grid
# resolution is set automatically)
NXI_CYLINDRICAL = 15

PMAX_KIN = 0.8  # maximum momentum (in mc) of hottail grid
PSEP_KIN = 0.05 # high grid-density region for p<psep (for mode SUPERTHERMAL, psep=0 effectively)

# Available simulation modes
MODE_FLUID        = 1
MODE_SUPERTHERMAL = 2
MODE_KINETIC      = 3
MODE_ISOTROPIC    = 4
# The following is a "fully kinetic" mode, with the temperature initiated at
# T~0, as in the superthermal mode. This mode is NOT part of the DREAM paper
# and remains mostly untested.
MODE_SUPERTHERMAL_KINETIC = 5

INCLUDE_FLUID_HOTTAIL = True    # Include fluid generation rates for hot-tail (in MODE_FLUID)
KINETIC_IONIZATION = True       # Kinetic ionization in non-fluid mode
SAUTER_CONDUCTIVITY = True      # Use the Sauter conductivity formula
SOLVER_RELTOL = 2e-6            # Default relative tolerance for Newton solver
LINEAR_SOLVER = Solver.LINEAR_SOLVER_MKL    # Linear solver to use
#LINEAR_SOLVER = Solver.LINEAR_SOLVER_LU

EXTENSION = ''

"""
#########################################
# HELPER ROUTINES
#########################################
"""

def setname(mode, scenario, phase, toroidal=True, io='settings'):
    return outname(mode=mode, scenario=scenario, phase=phase, toroidal=toroidal, io=io)


def outname(mode, scenario, phase, toroidal=True, io='output'):
    """
    Returns the appropriate output name to use.
    """
    base = '{}_{}_scen{}_{}_{}'.format(PATTERN, modename(mode), scenario, io, phase)

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
    elif i == MODE_SUPERTHERMAL_KINETIC: return 'hybrid'
    else: return ''


"""
#########################################
# SIMULATION ROUTINES
#########################################
"""

def getBaseline(mode=MODE_KINETIC, scenario=0, prescribedJ=False, toroidal=True, runInit=True, verboseInit=False):
    """
    Build baseline scenario.

    1. Run calculation of conductivity (only requires one step, dummy E). The
       conductivity is then used to determine the initial E profile required
       to achieve the desired initial current density profile.
    2. Run a simulation to set the desired current density profile by
       prescribing the electric field profile calculated in step 1.

    :param int mode:         Simulation mode to run in (aka model to use).
    :param int scenario:     Index of scenario to run (one scenario = one set of
                             disruption parameters).
    :param bool prescribedJ: If ``True``, sets the current density profile
                             returned by the function
                             ``Tokamak.getCurrentDensity()``. Otherwise assumes
                             a steady-state, fully ohmic plasma (i.e. constant
                             applied loop voltage)
    :param bool toroidal:    If ``True``, sets up a simulation in
                             toroidal/tokamak geometry. Otherwise, uses
                             cylindrical geometry.
    :param bool runInit:     If ``True``, runs the initialization simulation
                             that generates the desired current density profile.
                             This can be skipped if such an initialization
                             simulation has previously been run for the same
                             parameters.
    :param bool verboseInit: If ``True``, shows full output from DREAM when
                             running the initialization simulation.

    :return: Resulting DREAMSettings object.
    """
    lsetname = lambda phase : setname(mode=mode, scenario=scenario, phase=phase, toroidal=toroidal)
    loutname = lambda phase : outname(mode=mode, scenario=scenario, phase=phase, toroidal=toroidal)

    ds = DREAMSettings()

    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS #_COLLISIONAL

    # Set radial grid
    Tokamak.setMagneticField(ds, nr=NR, visualize=False, rGridNonuniformity=1, toroidal=toroidal)

    # Set dummy electric field
    ds.eqsys.E_field.setPrescribedData(1e-4)

    # Set temperature profile
    rT, T0 = Tokamak.getInitialTemperature()
    ds.eqsys.T_cold.setPrescribedData(T0, radius=rT)

    # Background ion density
    ds.eqsys.n_i.addIon(name='D', Z=1, Z0=1, iontype=Ions.IONS_DYNAMIC, T=T0, n=Tokamak.ne0*np.ones((rT.size,)), r=rT)

    # Background free electron density from ions
    nfree, rn0 = ds.eqsys.n_i.getFreeElectronDensity()

    # Use Sauter formula for conductivity?
    if SAUTER_CONDUCTIVITY:
        ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)
    
    # Disable kinetic grids during conductivity simulation
    ds.runawaygrid.setEnabled(False)
    ds.hottailgrid.setEnabled(False)
    
    # First calculation only requires one infinitesimal time step
    # (because we only need DREAM to evaluate the conductivity formula
    # with the prescribed temperature and density)
    ds.timestep.setTmax(1e-11)
    ds.timestep.setNt(1)

    # Calculate and store every fluid and scalar 'OtherQuantity' to extract as
    # much data as possible from the simulation.
    ds.other.include('fluid', 'scalar')

    # Calculate conductivity and radial grid
    do = DREAM.runiface(ds, quiet=True)


    ###############################################
    # PART 2
    # Generate initial current
    ####

    # Set kinetic grid?
    if mode != MODE_FLUID:
        ds.hottailgrid.setEnabled(True)
        # We initialize all distribution functions using a full linearized
        # collision operator in order to achieve the desired initial current
        # density.
        ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

        # Set non-uniform pitch grid
        if mode == MODE_ISOTROPIC:
            ds.hottailgrid.setNxi(1)
        elif not toroidal:
            ds.hottailgrid.setNxi(NXI_CYLINDRICAL)
        else:
            # Old way of settings trapped/passing xi grid
            #xi0Trapped = do.grid.xi0TrappedBoundary[:]
            #ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(xi0Trapped, dxiMax = 2/NR)

            # New way of setting trapped/passing xi grid
            ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax = 2/NR)

        # Set p grid
        ds.hottailgrid.setPmax(PMAX_KIN)
        if mode == MODE_KINETIC or mode == MODE_SUPERTHERMAL_KINETIC:
            ds.hottailgrid.setNp(NP_KIN1+NP_KIN2)
            ds.hottailgrid.setBiuniformGrid(psep=PSEP_KIN, npsep=NP_KIN1)
        elif mode == MODE_SUPERTHERMAL or mode == MODE_ISOTROPIC:
            ds.hottailgrid.setNp(NP_KIN2)

        mod = 1
        if mode == MODE_SUPERTHERMAL or mode == MODE_ISOTROPIC:
            # Scale density of initial Maxwellian so that there is always
            # a trace cold population (i.e. ncold > 0)
            mod = 0.9999
            # Corrected conductivity is only needed when the Maxwellian is
            # contained on the momentum grid (i.e. in MODE_KINETIC and
            # MODE_SUPERTHERMAL_KINETIC)
            ds.eqsys.j_ohm.setCorrectedConductivity(False)

        ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=mod*nfree, rT0=rT, T0=T0)
        # All electrons between pThreshold and pMax are counted as
        # hot electrons
        ds.eqsys.f_hot.setHotRegionThreshold(5)

        # Boundary condition on f at p = pMax (assume f(p>pMax) = 0)
        ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
        # Enable flux limiters
        ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)
        # Do not include the Jacobian elements for d f_hot / d n_i, i.e.
        # derivatives with respect to the ion densities (which would take
        # up *significant* space in the matrix)
        ds.eqsys.f_hot.enableIonJacobian(False)

    # Determine electric field needed for desired current
    # density profile
    if prescribedJ:
        # Prescribed current density profile
        rj, j = Tokamak.getCurrentDensity(r=do.grid.r[:])
        j  /= Tokamak.j0
        j0 = Tokamak.Ip * 2.0*np.pi / do.grid.integrate(j)
        print('Central plasma current density: {:.2f} MA/m^2'.format(j0/1e6))
        E0 = j0*j / do.other.fluid.conductivity[-1,:] * np.ones((1,rj.size))
        ds.eqsys.E_field.setPrescribedData(E0, radius=rj, times=[0])
    else:
        # Assume steady-state, fully ohmic plasma
        sigma = do.other.fluid.conductivity[-1,:]
        R0 = Tokamak.R0 + do.grid.r[-1]     # major radius coordinate of last r point
        R  = Tokamak.R0 + do.grid.r[:]
        E0 = 2*np.pi * Tokamak.Ip / (R0 * do.grid.integrate(sigma / R))
        E  = R0*E0 / R * np.ones((1,do.grid.r.size))

        ds.eqsys.E_field.setPrescribedData(E, radius=do.grid.r, times=[0])

    # Done with file, so close it...
    do.close()

    # Run a long initial simulation to obtain a steady-state
    # solution to E (desired if prescribed == False)
    ds.timestep.setTmax(1) 
    ds.timestep.setNt(3)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(LINEAR_SOLVER)
    ds.solver.setMaxIterations(100)
    ds.solver.tolerance.set(reltol=1e-7)
    
    if mode != MODE_FLUID:
        # Adjust absolute tolerances for quantities which can have
        # difficulties converging early during the simulation due
        # to being negligibly small
        ds.solver.tolerance.set('n_hot', abstol=1e5)
        ds.solver.tolerance.set('j_hot', abstol=1)

    ds.save(lsetname('init'))

    # Generate current
    INITFILE = loutname('init')
    if runInit:
        print('First current simulation')
        
        ds.solver.setVerbose(verboseInit)
        do = DREAM.runiface(ds, INITFILE, quiet=(not verboseInit))

        # Rescale to obtain exact total current
        E = do.eqsys.E_field[-1,:] * Tokamak.Ip / do.eqsys.I_p[-1]
        ds.eqsys.E_field.setPrescribedData(E, radius=do.grid.r, times=[0])

        # Close output
        do.close()

        print('Second current simulation')
        do = DREAM.runiface(ds, INITFILE, quiet=(not verboseInit))

    ##################################################
    # PART 3
    # Set up final DREAMSettings object to return
    ####
    # Copy settings object
    ds1 = DREAMSettings(ds)
    ignorelist = ['n_i', 'N_i', 'W_i']
    # Enable runaway generation
    ds1.eqsys.n_re.setEceff(RunawayElectrons.COLLQTY_ECEFF_MODE_FULL)
    if mode == MODE_FLUID:
        # Enable desired runaway generation rates
        ds1.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_FLUID_HESSLOW)
        ds1.eqsys.n_re.setDreicer(RunawayElectrons.DREICER_RATE_NEURAL_NETWORK)
        if INCLUDE_FLUID_HOTTAIL:
            ds1.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree, rT0=rT, T0=T0)
            ds1.eqsys.n_re.setHottail(RunawayElectrons.HOTTAIL_MODE_ANALYTIC_ALT_PC)

        #ds1.eqsys.f_re.enableAnalyticalDistribution()
    else:
        # Use fluid avalanche for isotropic...
        if mode == MODE_ISOTROPIC:
            ds1.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_FLUID_HESSLOW)
        # ...and kinetic avalanche for superthermal and kinetic...
        else:
            ds1.eqsys.n_re.setAvalanche(RunawayElectrons.AVALANCHE_MODE_KINETIC, pCutAvalanche = 0.01)

        # Fluid Dreicer must be disabled whenever a
        # distribution function is evolved
        ds1.eqsys.n_re.setDreicer(RunawayElectrons.DREICER_RATE_DISABLED)
        # Enable flux limiters
        ds1.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)

        # Account for fast electron impact ionization
        # (i.e. use f_hot to evaluate ionization rates)
        if KINETIC_IONIZATION:
            ds1.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)
        
        if mode == MODE_SUPERTHERMAL or mode == MODE_ISOTROPIC:
            # Take the superthermal limit of the linearized collision operator
            ds1.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
        
        if mode == MODE_SUPERTHERMAL or mode == MODE_ISOTROPIC or mode == MODE_SUPERTHERMAL_KINETIC:
            ds1.eqsys.T_cold.setInitialProfile(temperature=1)
            # Do not load initial values for these quantities
            # from the previous simulation because...
            ignorelist.append('T_cold')     # ...we want to start this at T=1
            # ...these are defined differently in MODE_SUPERTHERMAL and
            # MODE_ISOTROPIC and should therefore be recalculated when
            # launching the next simulation
            ignorelist.append('n_hot')
            ignorelist.append('n_cold')

        if mode == MODE_KINETIC or mode == MODE_SUPERTHERMAL_KINETIC:
            # Include a particle source/sink which handles creation/destruction
            # of electrons due to ionization/recombination
            ds1.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_EXPLICIT, shape=FHot.PARTICLE_SOURCE_SHAPE_DELTA)

    
    # Set self-consistent electric field and temperature
    ds1.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
    # Perfectly conducting wall
    ds1.eqsys.E_field.setBoundaryCondition(
        ElectricField.BC_TYPE_PRESCRIBED,
        inverse_wall_time=0, V_loop_wall_R0=0)
    ds1.eqsys.T_cold.setType(Temperature.TYPE_SELFCONSISTENT)
  
    # Set relative and absolute tolerances for hot and runaway
    # electron quantities which may sometimes be negligible and
    # therefore have difficulties converging
    MAXIMUM_IGNORABLE_ELECTRON_DENSITY = 1e5 #m-3
    ds1.solver.tolerance.set(reltol=SOLVER_RELTOL)
    ds1.solver.tolerance.set(unknown='n_re', reltol=SOLVER_RELTOL, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
    ds1.solver.tolerance.set(unknown='j_re', reltol=SOLVER_RELTOL, abstol=1e-10*MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
    if mode != MODE_FLUID:
        ds1.solver.tolerance.set(unknown='f_hot', reltol=SOLVER_RELTOL, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
        ds1.solver.tolerance.set(unknown='n_hot', reltol=SOLVER_RELTOL, abstol=MAXIMUM_IGNORABLE_ELECTRON_DENSITY)
        ds1.solver.tolerance.set(unknown='j_hot', reltol=SOLVER_RELTOL, abstol=1e-10*MAXIMUM_IGNORABLE_ELECTRON_DENSITY)

    if mode != MODE_FLUID:
        ds1.solver.tolerance.set(unknown='f_hot', reltol=SOLVER_RELTOL, abstol=1e5)
        ds1.solver.tolerance.set(unknown='n_hot', reltol=SOLVER_RELTOL, abstol=1e5)
        ds1.solver.tolerance.set(unknown='j_hot', reltol=SOLVER_RELTOL, abstol=1e-5)

    # Include information about time spent in different
    # parts of the code...
    ds1.output.setTiming(True, True)

    # Start from the state obtained in the init simulation
    ds1.fromOutput(INITFILE, ignore=ignorelist)
    return ds1


def simulate(ds1, mode, scenario, t_ioniz, nt_ioniz, t_TQ, nt_TQ, fracDinj, fracAr, Drr=0, runIoniz=True, runTQ=True, toroidal=True, verboseIoniz=False, verboseTQ=False):
    """
    Run the disruption simulation (given a baseline simulation). This function
    will inject a mixture of neutral deuterium and argon gas and run two phases.
    In the first phase, the ions are allowed to relax close to an equilibrium
    charge state distribution (occuring over a very short, microsecond, time
    scale). In the second phase the bulk of the thermal and current quenches
    occur.

    :param ds1:          Baseline simulation settings.
    :param mode:         Simulation mode/model (e.g. MODE_FLUID, MODE_ISOTROPIC,
                         MODE_SUPERTHERMAL, MODE_KINETIC or
                         MODE_SUPERTHERMAL_KINETIC).
    :param scenario:     ID of the scenario (i.e. set of plasma parameters)
                         to simulate.
    :param t_ioniz:      Simulation time of the ionization simulation.
    :param nt_ioniz:     Number of time steps to take in the ionization
                         simulation.
    :param t_TQ:         Simulation time of the thermal quench simulation.
    :param nt_TQ:        Number of time steps to take in thermal quench
                         simulation.
    :param fracDinj:     Amount of neutral deuterium to inject (fraction of
                         initial electron density).
    :param fracAr:       Amount of neutral argon to inject (fraction of initial
                         electron density).
    :param Drr:          Electron heat diffusion coefficient.
    :param runIoniz:     Run the ionization simulation (else, load from
                         previously executed ionization simulation).
    :param runTQ:        Run the thermal quench simulation (else, load from
                         previously executed thermal quench simulation).
    :param toroidal:     Use toroidal/tokamak geometry.
    :param verboseIoniz: Show output from DREAM during ionization simulation.
    :param verboseTQ:    Show output from DREAM during thermal quench simulation.
    """
    lsetname = lambda phase : setname(mode=mode, scenario=scenario, phase=phase, toroidal=toroidal)
    loutname = lambda phase : outname(mode=mode, scenario=scenario, phase=phase, toroidal=toroidal)

    ###################################
    # 1. IONIZATION PART OF TQ (~1 µs)
    ####
    # Add injected impurities
    ds1.eqsys.n_i.addIon('Dinj',  Z=1,  iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=fracDinj*Tokamak.ne0, T=1)
    ds1.eqsys.n_i.addIon('Ar',    Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=fracAr*Tokamak.ne0, T=1)
    
    # Prescribe heat diffusion?
    if Drr > 0:
        ds1.eqsys.T_cold.transport.prescribeDiffusion(drr=Drr)

    ds1.timestep.setTmax(t_ioniz)
    ds1.timestep.setNt(nt_ioniz)

    ds1.solver.setVerbose(verboseIoniz)
    ds1_outname = loutname('1')
    ds1.output.setFilename(ds1_outname)
    ds1.save(lsetname('1'))
    if runIoniz:
        DREAM.runiface(ds1,ds1_outname)

    ##################################
    # 2. THERMAL AND CURRENT QUENCH
    ####

    ds2 = DREAMSettings(ds1)
    ds2.fromOutput(ds1_outname)
    ds2.clearIgnore()

    # New time step and tMax
    ds2.timestep.setTmax(t_TQ - t_ioniz)
    ds2.timestep.setNt(nt_TQ)
    #ds2.timestep.setNumberOfSaveSteps(min(nt_TQ,1000))

    # In case the PARDISO (MKL) solver fails, we can switch
    # to the slightly more robust (but slower) built-in
    # PETSc linear solver.
    if LINEAR_SOLVER != Solver.LINEAR_SOLVER_MKL:
        ds2.solver.setBackupSolver(Solver.LINEAR_SOLVER_LU)

    # Set flux limiters
    if (scenario == 3 or scenario == 2) and mode == MODE_SUPERTHERMAL:
        ds2.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_UPWIND, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)
    elif scenario == 0 and mode == MODE_KINETIC:
        #ds2.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)
        ds2.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_UPWIND, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)

    ds2.solver.setVerbose(verboseTQ)
    ds2_outname = loutname('2')
    ds2.output.setFilename(ds2_outname)
    ds2.save(lsetname('2'))
    if runTQ:
        DREAM.runiface(ds2,ds2_outname)


def main():
    """
    Program entry point.
    """
    global EXTENSION

    # List of simulation scenarios. For the DREAM paper, scenarios
    # 1 and 4 were used (described in sections 6.2 and 6.3 respectively)
    scenarios = [
        # 0
        {'t_ioniz': 5e-6, 't_TQ': 4e-3, 'nt_ioniz': 400, 'nt_TQ': 2000, 'fracDinj': 10, 'fracAr': 0.6, 'Drr': 4000},
        # 1 ("full-conversion scenario")
        {'t_ioniz': 2e-5, 't_TQ': 1e-3,  'nt_ioniz': 200, 'nt_TQ': 1500, 'fracDinj': 1, 'fracAr': 1, 'Drr': 4000},
        # 2
        {'t_ioniz': 5e-6, 't_TQ': 9e-3, 'nt_ioniz': 400, 'nt_TQ': 1600, 'fracDinj': 10, 'fracAr': 0.15, 'Drr': 4000},
        # 3
        {'t_ioniz': 5e-6, 't_TQ': 4e-3, 'nt_ioniz': 1600, 'nt_TQ': 1500, 'fracDinj': 10, 'fracAr': 0.3, 'Drr': 4000},
        # 4 ("slow disruption scenario")
        {'t_ioniz': 1e-5, 't_TQ': 5e-3, 'nt_ioniz': 400, 'nt_TQ': 1600, 'fracDinj': 2, 'fracAr': 0.2, 'Drr': 1000}
    ]

    parser = argparse.ArgumentParser(description="Run a DREAM disruption simulation")
    parser.add_argument('-c', '--cylindrical', help="Run in cylindrical geometry", dest="cylindrical", action="store_true")
    parser.add_argument('-e', '--extension', help="Append the given string to the end of all file names", dest="extension", action="store", type=str)
    parser.add_argument('--fluid', help="Run simulation in pure fluid mode", dest='mode', action='store_const', const=MODE_FLUID)
    parser.add_argument('--kinetic', help="Run simulation in fully kinetic mode", dest='mode', action='store_const', const=MODE_KINETIC)
    parser.add_argument('--superthermal', help="Run simulation in superthermal mode", dest='mode', action='store_const', const=MODE_SUPERTHERMAL)
    parser.add_argument('--isotropic', help="Run simulation in isotropic mode", dest='mode', action='store_const', const=MODE_ISOTROPIC)
    parser.add_argument('--hybrid', help="Run simulation in hybrid superthermal-kinetic mode", dest='mode', action='store_const', const=MODE_SUPERTHERMAL_KINETIC)
    parser.add_argument('--steady-state-E', help="Assume that a steady state electric field is applied for initial current density", dest="steadyE", action='store_true')
    parser.add_argument('-r', '--run-from', help="Determines which simulation to start running from (0 = current, 1 = ionization, 2 = CQ", dest="runfrom", action='store', type=int)
    parser.add_argument('-s', '--scenario', help="Index of scenario to run 0-{}".format(len(scenarios)-1), dest='scenario', action='store', type=int)
    parser.add_argument('-v', '--verbose', help="Select which simulations to make verbose (0 = current, 1 = ionization, 2 = CQ", dest="verbose", nargs='*', type=int)

    parser.set_defaults(mode=MODE_FLUID, runfrom=0, scenario=0, verbose=[], extension='')

    # Parse arguments
    args = parser.parse_args()

    if args.extension:
        EXTENSION = '_' + args.extension

    # Select scenario (i.e. get scenario specific parameters)
    scenario = scenarios[args.scenario]
    # Set up baseline scenario
    runInit = (args.runfrom <= 0)
    verboseInit = (0 in args.verbose)
    ds = getBaseline(mode=args.mode, scenario=args.scenario, prescribedJ=not args.steadyE, toroidal=not args.cylindrical, runInit=runInit, verboseInit=verboseInit)

    # Run the simulation
    simulate(ds,
        mode=args.mode, scenario=args.scenario,
        runIoniz=(args.runfrom <= 1),
        runTQ=(args.runfrom<=2),
        toroidal=not args.cylindrical,
        verboseIoniz=(1 in args.verbose),
        verboseTQ=(2 in args.verbose),
        **scenario)

    return 0


if __name__ == '__main__':
    sys.exit(main())


