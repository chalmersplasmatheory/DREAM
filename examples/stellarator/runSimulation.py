import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append('../../py/')


from DREAM import DREAMSettings, DREAMOutput, DREAMException, DREAMQuitException, runiface

import DREAM.Settings.RadialGrid as RadialGrid
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
import DREAM.Settings.Equations.BootstrapCurrent as BootstrapCurrent

import stellarator

nr=101

activated=True

# First, download the configurations from Stefan Buller's website:
# https://sb0095.mycpanel.princeton.edu/QA/notes.html
eqfile = "EQs/QA_A5.nc"
if os.path.exists(f'{eqfile[:-3]}_settingsData.h5'):
    load = True
    print('yes')
else: 
    load = False
    print('no')


# Temperature evolution parameters 
tau = 0.5e-3 # decay time scale [ms]
Tend = 1     # Post-TQ temperature [eV]

# Current settings
bootstrap = True
Ip0 = 4e6

# Time parameters
tmax = tau*Tend*30 
nt = 10000

ds = DREAMSettings()
stellarator.setStellaratorRadialGrid(ds, eqfile, nr, load=load, fb=1.1)

if activated:
    rni, ni = stellarator.getInitialIonDensityProfile(1e20)
    ds.eqsys.n_i.addIon('D', n=ni, Z=1, Z0=1, r=rni, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
    ds.eqsys.n_i.addIon('T', n=ni, Z=1, Z0=1, r=rni, iontype=Ions.IONS_DYNAMIC, opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE, tritium=True)
else:
    rni, ni = stellarator.getInitialElectronDensityProfile(1e20)
    ds.eqsys.n_i.addIon('D', n=ni, Z=1, Z0=1, r=rni, iontype=Ions.IONS_DYNAMIC,opacity_mode=Ions.ION_OPACITY_MODE_GROUND_STATE_OPAQUE)
    
# Set collision settings
ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL
ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

# Use Sauter formula for conductivity
ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)

# Kinetic simulations not yet possible in stellarator mode
ds.runawaygrid.setEnabled(False)
ds.hottailgrid.setEnabled(False)

# Set initial current density 
rj, j = stellarator.getInitialCurrentDensityProfile(ds, bootstrap)
if bootstrap:
    ds.eqsys.j_bs.setMode(BootstrapCurrent.BOOTSTRAP_MODE_STELLARATOR)
    ds.eqsys.j_bs.setInitMode(BootstrapCurrent.BOOTSTRAP_INIT_MODE_TOTAL)
    ds.eqsys.j_ohm.setInitialProfile(j, radius=rj)
else:
    ds.eqsys.j_ohm.setInitialProfile(j, radius=rj, Ip0=Ip0)

# Set self consistent electric field evolution
ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_SELFCONSISTENT, inverse_wall_time=stellarator.tau_inv, R0=stellarator.R0)

# Set prescribed temperature evolution
rTe, tTe, Te = stellarator.getElectronTemperatureProfileEvolution(tau=tau, Tend=Tend)
ds.eqsys.T_cold.setType(Temperature.TYPE_PRESCRIBED)
ds.eqsys.T_cold.setPrescribedData(Te, radius=rTe, times=tTe)

nfree, rnfree = ds.eqsys.n_i.getFreeElectronDensity()
ds.eqsys.f_hot.setInitialProfiles(rn0=rnfree, n0=nfree, rT0=rTe, T0=Te[0,:])

# Set runaway electron settings
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)
ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_ANALYTIC_ALT_PC)

if activated:
    ds.eqsys.n_re.setCompton(Runaways.COMPTON_RATE_ITER_DMS)
    ds.eqsys.n_re.setTritium(Runaways.TRITIUM_MODE_FLUID)
else:
    ds.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_NEGLECT)
    ds.eqsys.n_re.setTritium(False)

ds.other.include(['fluid', 'scalar'])

# Set solver settings
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
ds.solver.setMaxIterations(maxiter=500)
ds.solver.tolerance.set(reltol=2e-6)
ds.solver.tolerance.set(unknown='n_re', reltol=2e-6, abstol=1e5)
ds.solver.tolerance.set(unknown='j_re', reltol=2e-6, abstol=1e-5)

ds.timestep.setNt(nt)
ds.timestep.setTmax(tmax)

ds.save('settings.h5')

do = runiface(ds, 'output.h5')
