#!/usr/bin/env python3
#
# This example demonstrates the use of consistent radial transport
# on f_re and n_re.
#
# Run as
#
#   $ ./generate.py
#   $ ../../build/iface/dreami dream_settings.h5
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport

ds = DREAMSettings()

# Physical parameters
E = 0.5     # Electric field strength (V/m)
n = 5e19    # Electron density (m^-3)
T = 1e3     # Temperature (eV)

# Grid parameters
pMax = 50     # maximum momentum in units of m_e*c
Np   = 50     # number of momentum grid points
Nxi  = 24     # number of pitch grid points
tMax = 1      # simulation time in seconds
Nt   = 50     # number of time steps
Nr   = 8      # number of radial grid points

dBOverB = 1e-3  # Magnetic perturbation strength
R0 = 1.6      # Tokamak major radius

# Set E_field
#ds.eqsys.E_field.setPrescribedData(E)
ds.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
ds.eqsys.E_field.setBoundaryCondition(bctype=Efield.BC_TYPE_PRESCRIBED, V_loop_wall_R0=E*2*np.pi, R0=R0)

# Set temperature
ds.eqsys.T_cold.setPrescribedData(T)

# Set ions
ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

# Disable avalanche generation
ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

# Disable hot-tail grid
ds.hottailgrid.setEnabled(False)
# Runaway grid
ds.runawaygrid.setNxi(Nxi)
ds.runawaygrid.setNp(Np)
ds.runawaygrid.setPmax(pMax)

# Set initial hot electron Maxwellian
#ds.eqsys.f_re.setInitialValue(0)

# Set up radial grid
ds.radialgrid.setB0(3)
ds.radialgrid.setMinorRadius(0.5)
ds.radialgrid.setWallRadius(0.5)
ds.radialgrid.setNr(Nr)

# Set Rechester-Rosenbluth transport
ds.eqsys.f_re.transport.setMagneticPerturbation(dBOverB)
ds.eqsys.f_re.transport.setBoundaryCondition(Transport.BC_F_0)

# Set solver type
#ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
ds.solver.setType(Solver.NONLINEAR)
ds.solver.setVerbose(False)
ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
#ds.solver.tolerance.set(reltol=1e-4)

# include otherquantities to save to output
ds.other.include('fluid')

# Set time stepper
ds.timestep.setTmax(tMax)
ds.timestep.setNt(Nt)

ds.output.setTiming(stdout=True, file=True)
ds.output.setFilename('output.h5')

# Save settings to HDF5 file
ds.save('dream_settings.h5')

