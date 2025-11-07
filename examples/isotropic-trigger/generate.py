#!/usr/bin/env python3
#
# Simple example which tests the implementation of the triggered
# isotropic model.
#

import matplotlib.pyplot as plt
import numpy as np
from DREAM import DREAMSettings, runiface
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.EquationTrigger as Trigger
import DREAM.Settings.Solver as Solver
import sys


def setup():
    ds = DREAMSettings()

    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS
    ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_SUPERTHERMAL

    ds.radialgrid.setMinorRadius(1.0)
    ds.radialgrid.setWallRadius(1.2)
    ds.radialgrid.setNr(20)
    ds.radialgrid.setB0(5.6)

    # Radial grid for input parameters
    r0 = np.linspace(0, ds.radialgrid.a, ds.radialgrid.nr)
    rho = r0 / ds.radialgrid.a

    # Electric field
    ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(
        EField.BC_TYPE_PRESCRIBED,
        inverse_wall_time=0, V_loop_wall_R0=0
    )
    j0 = (1 - (1-0.001**(1/0.41))*rho**2)**0.41
    ds.eqsys.j_ohm.setInitialProfile(j0, r0, Ip0=15e6)

    ds.hottailgrid.setEnabled(True)
    ds.runawaygrid.setEnabled(False)

    ds.hottailgrid.setPmax(0.8)
    ds.hottailgrid.setNp(80)
    ds.hottailgrid.setNxi(1)

    # Set initial temperature
    T0 = 20e3 * (1 - 0.99*rho**2)
    n0 = 1e20 * np.ones(r0.shape)
    ds.eqsys.T_hot.setInitialProfile(T0, r0)

    # Add ion species
    ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n0, r=r0)
    ds.eqsys.n_i.addIon('Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=0)

    # Set initial properties of f_hot
    mod = 0.9999
    ds.eqsys.f_hot.setInitialProfiles(mod*n0, T0, rn0=r0, rT0=r0)
    ds.eqsys.f_hot.trigger.equation.enableIonJacobian(False)
    ds.eqsys.f_hot.trigger.equation.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)

    # Enable trigger
    #ds.eqsys.f_hot.enableIsotropicTrigger(Trigger.TYPE_COLD_ELECTRON_RISE)
    ds.eqsys.f_hot.enableIsotropicTrigger(Trigger.TYPE_TIME, trigger_time=1)

    # Solver settings
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
    ds.solver.setMaxIterations(100)

    ds.solver.tolerance.set(reltol=1e-6)
    ds.solver.tolerance.set('j_hot', abstol=1)

    # Time step settings
    ds.timestep.setTmax(1e-5)
    ds.timestep.setNt(100)

    return ds


def main():
    ds = setup()
    ds.save('settings.h5')

    runiface(ds, 'output.h5')

    return 0


if __name__ == '__main__':
    sys.exit(main())


