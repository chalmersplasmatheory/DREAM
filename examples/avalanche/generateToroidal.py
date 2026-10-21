"""Demonstration of AVALANCHE_MODE_BOLTZMANN_KNOCK_ON.

A more advanced example in toroidal magnetic field and with
both hot-tail and runaway grids modelled kinetically.

The resulting growth rate is not fully resolved;
increased pmax, Np, Nxi and number of time steps will affect the results.

Run simulation with:
- python generateToroidal.py
- dreami settingsToroidal.h5
Visualize the results with plotAdvanced.py (requires the generateAdvanced.py
case to be run as well).

Approximate simulation time: 1 minute
"""

import numpy as np
import scipy.constants
import sys
from pathlib import Path

from generateAdvanced import pLong, generate_pitch_grid
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
kappa = 1.5  # elonggation
GOverR0 = 5.0  # toroidal magnetic field on axis
rref = np.linspace(0, a, 20)
IpRef = 1e6  # reference plasma current which generates the poloidal magnetic field (assumed uniform profile)
psiref = -mu0 * IpRef * (1 - (rref / a) ** 2) * a

rDelta = np.linspace(0, a, 20)
Delta = np.linspace(0, 0.1 * a, rDelta.size)
rdelta = np.linspace(0, a, 20)
delta = np.linspace(0, 0.2, rdelta.size)

#########################
# RESOLUTION PARAMETERS #
#########################

nTimeSteps = 150
tmax = 0.3

# Hot grid
Nxi = 30
npsep = 50
nxiexp = 1.5
psep = 1.2
pmax = 30
Np = 250


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
ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)
ds.eqsys.f_hot.setBoundaryCondition(FHot.BC_F_0)
ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_IMPLICIT)

ds.eqsys.n_re.setAvalanche(
    avalanche=Runaways.AVALANCHE_MODE_BOLTZMANN_KNOCK_ON,
    pCutAvalanche=0.1,
)
# ds.eqsys.n_re.setEceff(Eceff=Runaways.COLLQTY_ECEFF_MODE_SIMPLE)
ds.eqsys.n_re.setInitialProfile(
    density=0.01
)  # arbitrary initial value for n_re to seed the avalanche
ds.eqsys.f_hot.enableIonJacobian(False)
ds.runawaygrid.setEnabled(False)

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

phot_f = pLong[pLong<=pmax]
ds.hottailgrid.setCustomGrid(p_f=phot_f, xi_f=generate_pitch_grid(Nxi, nxiexp))

ds.timestep.setTmax(tmax)
ds.timestep.setNt(nTimeSteps)

ds.output.setFilename("outputToroidal.h5")
ds.output.setTiming(stdout=True)

ds.other.include("fluid", "hottail/C_boltz", "hottail/S_ava")

ds.save("settingsToroidal.h5")
