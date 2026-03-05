"""Demonstration of AVALANCHE_MODE_BOLTZMANN_KNOCK_ON.

A more advanced example in toroidal magnetic field and with
both hot-tail and runaway grids modelled kinetically.

The resulting growth rate is not fully resolved;
increased pmax, Np, Nxi and number of time steps will affect the results.

For stability reasons, the simulation is set up as a 3-stage rocket:
 - python generateAdvanced.py
 - dreami settingsAdvanced1.h5
 - dreami settingsAdvanced2.h5
 - dreami settingsAdvanced3.h5
Each simulation load the previous solution's results and build to
a quasi-steady state growth rate.

Visualize solutions with plotAdvanced.py

Approximate simulation time: 5 minutes
"""

import numpy as np
import scipy.constants
import sys
from pathlib import Path

try:
    import DREAM
except ModuleNotFoundError:
    sys.path.append(str(Path(__file__).resolve().parents[2] / "py"))
    import DREAM


from DREAM.DREAMSettings import DREAMSettings

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.IonSpecies as IonSpecies
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver


######################
# PHYSICAL CONSTANTS #
######################
c = scipy.constants.c  # speed of light
ec = scipy.constants.e  # elementary charge
eps0 = scipy.constants.epsilon_0  # vacuum permittivity
me = scipy.constants.m_e  # electron mass

################################
# SIMULATION PLASMA PARAMETERS #
################################
EOverEc = 10
T = 10
EOverEcTot = EOverEc
nD0 = 1e19
nD1 = 5e18
nAr = 5e18
nNe = 0

nTot = nD0 + nD1 + nAr * 18 + nNe * 10  # total (free plus bound) electron density
nFree = nD1  # free electron density

lnLambda = 14.9 - 0.5 * np.log(nFree / 1e20) + np.log(T / 1e3)
EcTot = nTot * lnLambda * (ec**3) / (4 * np.pi * (eps0**2) * me * (c**2))
E = EcTot * EOverEcTot

# ANALYTIC TOROIDAL MAGNETIC FIELD #
nr = 1
mu0 = scipy.constants.mu_0
R0 = 1  # major radius
a = 0.5  # minor radius
kappa = 1.5  # elongation
GOverR0 = 5.0  # toroidal magnetic field on axis
rref = np.linspace(0, a, 20)
IpRef = 1e6  # reference plasma current which generates the poloidal magnetic field (assumed uniform profile)
psiref = -mu0 * IpRef * (1 - (rref / a) ** 2) * a

rDelta = np.linspace(0, a, 20)
Delta = np.linspace(0, 0.1 * a, rDelta.size)
rdelta = np.linspace(0, a, 20)
delta = np.linspace(0, 0.2, rdelta.size)

# hardcoded TP boundary for this setup,
# needs to be adjusted if magnetic field changes.
# DREAM does not yet support non-uniform magnetic fields
# with automatic addition of boundary-layer cells.
xi0_TP_boundary = 0.62638244

def generate_pitch_grid(Nxi, nxiexp):
    # non-uniform pitch grid with geometric growth:
    #   xi=1:  dxi = 1 / Nxi^nxiexp
    #   xi=-1: dxi = nxiexp / Nxi
    # where we also add xi=0 as well as boundary-layer cells around xi = +/- xi0_TP_boundary.
    xi0_f = list(reversed(1 - 2 * np.linspace(0, 1, Nxi) ** nxiexp))
    if 0 not in xi0_f:
        xi0_f.append(0)
    for xi0T in [-xi0_TP_boundary, xi0_TP_boundary]:
        for k in [-1, 1, -3, 3]:
            xi0_f.append(xi0T + k * 1e-4)
    xi0_f = sorted(xi0_f)
    return xi0_f



def generate_momentum_grid(pMax, pCrit, relativeSize=70):
    """
    Creates a non-uniform momentum grid with quadratic grid spacing
        dp ~ p/relativeSize
    capped to dp(pCrit). Below 0.5 pCrit a coarser grid is used,
    assuming no interesting dynamics in that region.

    Grid boundaries are at p=0 and p=pMax
    """
    p_latest = 0
    ps = [p_latest]
    while  p_latest < pMax:
        if p_latest<0.5*pCrit:
            # coarse grid in low-speed region
            dp = 5*pCrit/relativeSize
        else:
            dp = max(pCrit/relativeSize, p_latest/relativeSize)
        p_latest += dp
        ps.append(p_latest)
    ps[-1] = pMax
    return np.array(ps, dtype=np.float64)

#########################
# RESOLUTION PARAMETERS #
#########################

nTimeSteps = 150
tmax = 0.3

# Hot grid
Nxi_hot = 30
npsep_hot = 50
nxiexp_hot = 1.2
psep_hot = 1.2
pmax_hot = 3
Np_hot = 100

# RE grid
Nxi_re = 50
nxiexp_re = 2
pmax_re = 240
Np_re = 500

###############
# DREAM SETUP #
###############
ds = DREAMSettings()

ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

ds.eqsys.E_field.setPrescribedData(E)

ds.eqsys.n_i.addIon(
    name="D_ionized", Z=1, n=nD0, iontype=IonSpecies.IONS_PRESCRIBED_FULLY_IONIZED
)
ds.eqsys.n_i.addIon(
    name="D_neutral", Z=1, n=nD1, iontype=IonSpecies.IONS_PRESCRIBED_NEUTRAL
)
ds.eqsys.n_i.addIon(name="Ar", Z=18, n=nAr, iontype=IonSpecies.IONS_PRESCRIBED_NEUTRAL)
ds.eqsys.n_i.addIon(name="Ne", Z=10, n=nNe, iontype=IonSpecies.IONS_PRESCRIBED_NEUTRAL)

ds.eqsys.T_cold.setPrescribedData(T)

# initialize f_hot to something small but smooth in order for the
# advection interpolation coefficients to converge but n_hot be
# negligible compared to n_re(t=0)
ds.eqsys.f_hot.setInitialProfiles(n0=1, T0=1e5)
ds.eqsys.f_hot.enableIonJacobian(False)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_IMPLICIT)

ds.runawaygrid.setEnabled(True)
#ds.eqsys.f_re.setInitialProfiles(n0=1, T0=1e5)
ds.eqsys.f_re.setAdvectionInterpolationMethod(FHot.AD_INTERP_TCDF)
ds.eqsys.f_re.setBoundaryCondition(FHot.BC_F_0)

ds.eqsys.n_re.setAvalanche(
    avalanche=Runaways.AVALANCHE_MODE_BOLTZMANN_KNOCK_ON,
    pCutAvalanche=0.1,
)
# ds.eqsys.n_re.setEceff(Eceff=Runaways.COLLQTY_ECEFF_MODE_SIMPLE)
# ds.eqsys.n_re.setInitialProfile(
#     density=1
# )  # arbitrary initial value for n_re to seed the avalanche

### Configure non-uniform magnetic field ###
ds.radialgrid.setShaping(
    psi=psiref,
    rpsi=rref,
    GOverR0=GOverR0,
    kappa=kappa,
    delta=delta,
    rdelta=rdelta,
    Delta=Delta,
    rDelta=rDelta,
)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(a)
ds.radialgrid.setMajorRadius(R0)
ds.radialgrid.setNr(nr)

ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
ds.solver.tolerance.set(reltol=1e-4)
ds.solver.setVerbose(True)

### Grids ###
# single big grid array p_f which we slice into a hot and re grid
pLong = generate_momentum_grid(pMax=pmax_re, pCrit=0.8, relativeSize=70)
idxRE = np.argwhere(pLong>=pmax_hot)[0][0]
phot_f = pLong[:idxRE+1]
pre_f = pLong[idxRE:] # shared endpoint with phot

ds.hottailgrid.setCustomGrid(p_f=phot_f, xi_f=generate_pitch_grid(Nxi_hot, nxiexp_hot))
ds.runawaygrid.setCustomGrid(p_f=pre_f, xi_f=generate_pitch_grid(Nxi_re, nxiexp_re))

ds.timestep.setTmax(tmax)
ds.timestep.setNt(nTimeSteps)

ds.output.setTiming(stdout=True)

ds.other.include("fluid", "hottail/C_boltz", "runaway/C_boltz")

if __name__=="__main__":
    ds_initial = DREAMSettings(ds, chain=False)

    # Step 1: create a runaway tail with UPWIND
    ds_initial.eqsys.f_hot.setAdvectionInterpolationMethod(FHot.AD_INTERP_UPWIND)
    ds_initial.eqsys.f_re.setAdvectionInterpolationMethod(FHot.AD_INTERP_UPWIND)
    ds_initial.timestep.setTmax(0.2)
    ds_initial.timestep.setNt(20)
    ds_initial.output.setFilename("outputAdvanced1.h5")
    ds_initial.save("settingsAdvanced1.h5")

    # Step 2: take a baby step to get a converged solution with flux limiter
    ds_restart1 = DREAMSettings(ds_initial, chain=True)
    ds_restart1.eqsys.f_hot.setAdvectionInterpolationMethod(FHot.AD_INTERP_TCDF)
    ds_restart1.eqsys.f_re.setAdvectionInterpolationMethod(FHot.AD_INTERP_TCDF)
    ds_restart1.timestep.setTmax(tmax/20)
    ds_restart1.timestep.setNt(nTimeSteps/5)
    ds_restart1.output.setFilename("outputAdvanced2.h5")
    ds_restart1.save("settingsAdvanced2.h5")

    # Step 3: big run to quasi-steady state
    ds_restart2 = DREAMSettings(ds_restart1, chain=True)
    ds_restart2.timestep.setTmax(tmax)
    ds_restart2.timestep.setNt(nTimeSteps)
    ds_restart2.output.setFilename("outputAdvanced.h5")
    ds_restart2.save("settingsAdvanced3.h5")

