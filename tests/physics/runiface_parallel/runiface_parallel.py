#!/usr/bin/env python3
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
import dreamtests

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
from DREAM.DREAMException import DREAMException
import DREAM

def createExperimentData(index):
    ds = DREAMSettings()
    #ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

    # Physical parameters
    E = 6       # Electric field strength (V/m)
    n = 5e19    # Electron density (m^-3)
    T = 100     # Temperature (eV)

    # Grid parameters
    pMax = 1    # maximum momentum in units of m_e*c
    Np   = 300  # number of momentum grid points
    Nxi  = 20   # number of pitch grid points
    tMax = 1e-3 # simulation time in seconds
    Nt   = 20   # number of time steps

    # Set E_field
    ds.eqsys.E_field.setPrescribedData(E)

    # Set temperature
    ds.eqsys.T_cold.setPrescribedData(T)

    # Set ions
    ds.eqsys.n_i.addIon(name='D', Z=1 + index, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

    # Disable avalanche generation
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

    # Hot-tail grid settings
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

    # Set initial hot electron Maxwellian
    ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)

    # Set boundary condition type at pMax
    #ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST) # extrapolate flux to boundary
    ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
    ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

    # Disable runaway grid
    ds.runawaygrid.setEnabled(False)

    # Set up radial grid
    ds.radialgrid.setB0(5)
    ds.radialgrid.setMinorRadius(0.22)
    ds.radialgrid.setWallRadius(0.22)
    ds.radialgrid.setNr(1)

    # Set solver type
    ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
    ds.solver.preconditioner.setEnabled(False)

    # include otherquantities to save to output
    ds.other.include('fluid','nu_s','nu_D')

    # Set time stepper
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(Nt)
    ds.output.setTiming(stdout=True, file=True)
    ds.output.setFilename(f'output{index}.h5')

    return ds

def run(args):
    number_of_tasks = 10
    ds = [createExperimentData(i) for i in range(number_of_tasks)]
    output = [None] * number_of_tasks
    try:
        result = DREAM.runiface_parallel(ds, output, quiet=True)
        dreamtests.print_ok("Parallel test run correctly.")
        return True
        # for r in result:
        #     print(r.eqsys.j_ohm[-1,:])
    except DREAMException as error:
        dreamtests.print_error(f"Parallel test failed: {error}")
        return False
