# ADAPTIVE TIME STEPPER TEST
#
# This test evolves a simple equation system a few steps in time using both the
# constant time stepper, as well as with the adaptive time stepper (but with a
# constant time step). If the results are not exactly the same, this indicates
# a problem with the adaptive time stepper.

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import sys

import dreamtests

import DREAM
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TimeStepper as TimeStepper
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions


def genSettings(adaptive=False):
    """
    Generate the baseline DREAMSettings object.
    """
    ds = DREAM.DREAMSettings()

    a    = 0.5
    B0   = 5
    E    = 10
    Nr   = 4
    Np   = 300
    Nt   = 5
    Nxi  = 5
    pMax = 0.2
    tMax = 1e-3
    T    = 50

    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL

    ds.radialgrid.setB0(B0)
    ds.radialgrid.setNr(Nr)
    ds.radialgrid.setMinorRadius(a)

    ds.timestep.setTmax(tMax)

    if adaptive:
        ds.timestep.setType(TimeStepper.TYPE_ADAPTIVE)
        ds.timestep.setDt(tMax / Nt)
        ds.timestep.setConstantStep(True)
    else:
        ds.timestep.setType(TimeStepper.TYPE_CONSTANT)
        ds.timestep.setNt(Nt)

    ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e20)
    ds.eqsys.E_field.setPrescribedData(E)
    ds.eqsys.T_cold.setPrescribedData(T)

    ds.hottailgrid.setEnabled(True)
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

    ds.runawaygrid.setEnabled(False)
    
    nfree_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree_initial, rT0=0, T0=T)
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_GMRES)
    ds.solver.setTolerance(reltol=0.01)

    return ds


def testrun():
    """
    Simple test run to make sure everything works.
    """
    ds = genSettings(True)
    ds.save('input.h5')
    do = DREAM.runiface(ds)

    do.eqsys.f_hot.plot(t=[0,2,4,5])
    plt.show()


def run(args):
    """
    Run the test.
    """
    #testrun()
    #return True


    QUIET = True
    TOLERANCE = 100*sys.float_info.epsilon

    output_const, output_adapt = None, None
    # Save output?
    if args['save']:
        output_const = 'output_const.h5'
        output_adapt = 'output_adapt.h5'

    # Run with constant time step
    ds_const = genSettings(False)
    do_const = DREAM.runiface(ds_const, output_const, quiet=QUIET)

    ds_adapt = genSettings(True)
    do_adapt = DREAM.runiface(ds_adapt, output_adapt, quiet=QUIET)

    # Compare results
    err_f = np.abs(do_const.eqsys.f_hot[-1,:] / do_adapt.eqsys.f_hot[-1,:] - 1)
    err_j = np.abs(do_const.eqsys.j_hot[-1,:] / do_adapt.eqsys.j_hot[-1,:] - 1)

    eps_f = np.amax(err_f)
    eps_j = np.amax(err_j)

    success = True
    if eps_f > TOLERANCE:
        dreamtests.print_error("f_hot is too different. eps_f = {:.8e}".format(eps_f))
        success = False
    if eps_j > TOLERANCE:
        dreamtests.print_error("j_hot is too different. eps_j = {:.8e}".format(eps_j))
        success = False

    if success:
        dreamtests.print_ok("Solution from adaptive time stepper agrees exactly with solution from the constant time stepper.")

    return success

