"""Demonstration of AVALANCHE_MODE_BOLTZMANN_KNOCK_ON.

A simple example in a uniform magnetic field with neutral impurities
and relatively weak electric field.

The resulting growth rate is not fully resolved;
increased pmax, Np, Nxi and number of time steps will affect the results.

Approximate simulation time: 10-30 seconds.
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

#########################
# RESOLUTION PARAMETERS #
#########################
Nxi = 30
npsep = 40
nxiexp = 1.5
nTimeSteps = 100
tmax = 0.2
psep = 0.7
pmax = 30
Np = 200

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

ds.radialgrid.setB0(1e-6)
ds.radialgrid.setMinorRadius(0.1)
ds.radialgrid.setWallRadius(0.1)
ds.radialgrid.setNr(1)

ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
ds.solver.tolerance.set(reltol=1e-4)
ds.solver.setVerbose(True)

# bi-uniform p grid
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pmax)
ds.hottailgrid.setBiuniformGrid(psep=psep, npsep=npsep)

# non-uniform pitch grid with geometric growth:
#   xi=1:  dxi = 1 / Nxi^nxiexp
#   xi=-1: dxi = nxiexp / Nxi

xi0_f = list(reversed(1 - 2 * np.linspace(0, 1, Nxi) ** nxiexp))
ds.hottailgrid.setCustomGrid(xi_f=xi0_f)

ds.timestep.setTmax(tmax)
ds.timestep.setNt(nTimeSteps)

ds.output.setFilename("outputCylindrical.h5")
ds.output.setTiming(stdout=True)

ds.other.include("fluid")
ds.other.include("hottail/C_boltz", "hottail/S_ava")

ds.save("settingsCylindrical.h5")
