#!/usr/bin/env python3
#SBATCH --job-name="Ware Pinch"
#SBATCH --output=logs/dream_%j.out
#SBATCH --error=logs/dream_%j.err
#SBATCH --time=40:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --chdir=/home/lowebl/DREAM/examples/warepinch
#
# This example shows how to set up a simple CODE-like runaway
# scenario in DREAM. The simulation uses a constant temperature,
# density and electric field, and generates a runaway current
# through the electric field acceleration, demonstrating Dreicer generation.
#
# Run as
#
#   $ ./basic.py
#   $ ../../build/iface/dreami dream_settings.h5
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Equations.HotElectronDistribution as f_hot

from DREAM import DREAMSettings, DREAMOutput, DREAMException, runiface

ds = DREAMSettings()
#ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

# Physical parameters
E = 6      # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 100     # Temperature (eV)
a = 0.22 #minor radius
b = 0.22 #minor wall radius
R0 = 0.68 #major radius
B0 = 5 #magnetic field at R0

# Grid parameters
pMax = 0.2   # maximum momentum in units of m_e*c
Np   = 100  # number of momentum grid points
Nxi  = 20   # number of pitch grid points
tMax = 4e-2 # simulation time in seconds
Nt   = 1600   # number of time steps
NR = 10 # Radial grid resolution

# Set E_field
ds.eqsys.E_field.setPrescribedData(E)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

# Hot-tail grid settings
ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax=1e-1) #sets trapping particles
ds.hottailgrid.setNp(Np)
ds.hottailgrid.setPmax(pMax)

# Set initial hot electron Maxwellian
ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)
ds.eqsys.f_hot.setParticleSource(f_hot.PARTICLE_SOURCE_IMPLICIT)

# Set boundary condition type at pMax
#ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST) # extrapolate flux to boundary
ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
ds.eqsys.f_hot.setWarePinchMode(DistFunc.WAREPINCH_MODE_NEGLECT)
ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_TCDF,ad_jac=DistFunc.AD_INTERP_JACOBIAN_UPWIND)
ds.eqsys.f_hot.transport.prescribeDiffusion(0.1)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
#ds.radialgrid.setB0(5)
#ds.radialgrid.setMinorRadius(0.22)
#ds.radialgrid.setWallRadius(0.22)
#ds.radialgrid.setNr(10) #radial resolution ten points

# Set up toroidal grid

# Radial grid for analytical magnetic field
r = np.linspace(0, a, NR)

# Elongation profile (fit)
kappa = np.polyval([0.13408883,-0.35950981,0.37803907,1.17567818], r)

# Triangularity profile
delta = 0.35*(r/a)

# Poloidal flux radial grid
psi_r = [0.,0.15729669,0.21956937,0.26821122,0.30894407,0.34456487,0.37669775,0.40612651,0.43344149,0.45903068,0.48318347,0.50611202,0.52798629,0.54893893,0.56907693,0.58848811,0.60724724,0.62541499,0.64304163,0.66016932,0.67683145,0.69305626,0.70886366,0.72427377,0.73930629,0.7539885,0.76834598,0.78240391,0.7961856,0.80971195,0.82300201,0.83607276,0.84893872,0.86161257,0.87410597,0.88642996,0.89859466,0.91060962,0.92248405,0.93422631,0.94584417,0.95734403,0.96873219,0.98001427,0.99119554,1.00228094,1.01327516,1.02418297,1.03500895,1.04575739,1.05643194,1.06703588,1.07757233,1.08804434,1.09845488,1.10880682,1.11910278,1.129345,1.13953555,1.14967642,1.15976976,1.16981777,1.17982258,1.18978613,1.19971007,1.20959591,1.21944516,1.22925925,1.23903956,1.24878743,1.25850418,1.26819105,1.27784926,1.28747997,1.29708427,1.30666306,1.31621725,1.32574772,1.33525552,1.34474181,1.3542077,1.36365421,1.37308209,1.38249198,1.39188451,1.40126024,1.41061957,1.41996291,1.42929063,1.43860331,1.44790164,1.45718634,1.46645808,1.47571723,1.48496403,1.49419869,1.50342147,1.51263288,1.52183344,1.53102371,1.54020416,1.54937504,1.55853657,1.56768898,1.57683252,1.58596753,1.59509438,1.60421339,1.61332461,1.62242774,1.63152244,1.64060842,1.64968555,1.65875389,1.6678135,1.67686424,1.6859037,1.69492821,1.7039341,1.71291765,1.72187424,1.73079888,1.73968684,1.74853379,1.75733735,1.76609582,1.77480767,1.78347367,1.79100379]
    # Poloidal flux
psi_p = [0.,0.02707751,0.05415503,0.08123254,0.10831006,0.13538757,0.16246508,0.1895426,0.21662011,0.24369763,0.27077514,0.29785265,0.32493017,0.35200768,0.3790852,0.40616271,0.43324022,0.46031774,0.48739525,0.51447276,0.54155028,0.56862779,0.59570531,0.62278282,0.64986033,0.67693785,0.70401536,0.73109288,0.75817039,0.7852479,0.81232542,0.83940293,0.86648045,0.89355796,0.92063547,0.94771299,0.9747905,1.00186802,1.02894553,1.05602304,1.08310056,1.11017807,1.13725559,1.1643331,1.19141061,1.21848813,1.24556564,1.27264315,1.29972067,1.32679818,1.3538757,1.38095321,1.40803072,1.43510824,1.46218575,1.48926327,1.51634078,1.54341829,1.57049581,1.59757332,1.62465084,1.65172835,1.67880586,1.70588338,1.73296089,1.76003841,1.78711592,1.81419343,1.84127095,1.86834846,1.89542598,1.92250349,1.949581,1.97665852,2.00373603,2.03081354,2.05789106,2.08496857,2.11204609,2.1391236,2.16620111,2.19327863,2.22035614,2.24743366,2.27451117,2.30158868,2.3286662,2.35574371,2.38282123,2.40989874,2.43697625,2.46405377,2.49113128,2.5182088,2.54528631,2.57236382,2.59944134,2.62651885,2.65359637,2.68067388,2.70775139,2.73482891,2.76190642,2.78898393,2.81606145,2.84313896,2.87021648,2.89729399,2.9243715,2.95144902,2.97852653,3.00560405,3.03268156,3.05975907,3.08683659,3.1139141,3.14099162,3.16806913,3.19514664,3.22222416,3.24930167,3.27637919,3.3034567,3.33053421,3.35761173,3.38468924,3.41176676,3.43884427,3.46592178]

ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
ds.radialgrid.setMajorRadius(R0)
ds.radialgrid.setMinorRadius(a)
ds.radialgrid.setWallRadius(b)
ds.radialgrid.setNr(NR)

ds.radialgrid.setShaping(psi=psi_p, rpsi=psi_r, GOverR0=B0, kappa=kappa, rkappa=r, delta=delta, rdelta=r)


# Set solver type
ds.solver.setType(Solver.NONLINEAR) # semi-implicit time stepping
ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
#ds.solver.setVerbose(True)
ds.solver.preconditioner.setEnabled(False)

# include otherquantities to save to output
#ds.other.include('hottail/Ar')
ds.other.include('fluid','nu_s','nu_D','hottail/Ar')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('ware_off_D=0.1.h5')

# Save settings to HDF5 file
ds.save('ware_off_D=0.1_settings.h5')

do = runiface(ds, 'ware_off_D=0.1.h5')
