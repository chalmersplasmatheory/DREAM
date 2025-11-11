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
import os
import sys


WITH_TRIGGER = True

MINOR_RADIUS = 1.0
NR = 20
TRIGGER_TIME = 1e-7


def set_geometry(ds):
    """
    Set the geometry of the simulation domain.
    """
    ds.radialgrid.setMinorRadius(MINOR_RADIUS)
    ds.radialgrid.setWallRadius(MINOR_RADIUS*1.2)
    ds.radialgrid.setNr(NR)
    ds.radialgrid.setB0(5.6)


def calculate_E(r, T, n, j, Ip):
    """
    Calculate the electric field profile, given temperature, density, and
    current density profiles. 
    """
    ds = DREAMSettings()

    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS
    ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_SUPERTHERMAL

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    set_geometry(ds)

    ds.eqsys.T_cold.setPrescribedData(T, radius=r)
    ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n, r=r)

    ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(
        EField.BC_TYPE_PRESCRIBED,
        inverse_wall_time=0, V_loop_wall_R0=0
    )

    ds.eqsys.j_ohm.setInitialProfile(j, r, Ip0=Ip)

    ds.timestep.setNt(1)
    ds.timestep.setTmax(1e-6)
    
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
    ds.output.setFilename('output_init.h5')

    do = runiface(ds, 'output_init.h5', quiet=True)

    return do.grid.r[:], do.eqsys.E_field[-1,:], ds


def setup():
    T0 = 20e3
    n0 = 1e20

    # Radial grid for input parameters
    r0 = np.linspace(0, MINOR_RADIUS, NR+1)
    rho = r0 / MINOR_RADIUS

    T = T0 * (1 - 0.99*rho**2)
    n = n0 * np.ones(r0.shape)

    j = (1-0.99*rho**2)
    rE, E, ds_init = calculate_E(r=r0, T=T, n=n, j=j, Ip=15e6)

    ds = DREAMSettings(ds_init)

    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda            = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode          = Collisions.PSTAR_MODE_COLLISIONLESS
    ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_SUPERTHERMAL

    set_geometry(ds)

    # Electric field
    ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(
        EField.BC_TYPE_PRESCRIBED,
        inverse_wall_time=0, V_loop_wall_R0=0
    )
    ds.eqsys.E_field.setInitialProfile(E, radius=rE)

    ds.hottailgrid.setEnabled(True)
    ds.runawaygrid.setEnabled(False)

    ds.hottailgrid.setPmax(1.4)
    ds.hottailgrid.setNp(80)
    ds.hottailgrid.setNxi(1)

    # Set initial temperature
    if WITH_TRIGGER:
        ds.eqsys.T_hot.setInitialProfile(T, r0)

    # Add ion species
    #ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=n, r=r0)
    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)
    if WITH_TRIGGER:
        ds.eqsys.n_i.addIon('Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=3e20)
    else:
        ds.eqsys.n_i.addIon('Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=3e20)

    # Set initial properties of f_hot
    mod = 0.9999
    ds.eqsys.f_hot.setInitialProfiles(mod*n, T, rn0=r0, rT0=r0)
    if WITH_TRIGGER:
        ds.eqsys.f_hot.trigger.equation.enableIonJacobian(False)
        ds.eqsys.f_hot.trigger.equation.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)
    else:
        ds.eqsys.f_hot.enableIonJacobian(False)
        ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF, ad_jac=FHot.AD_INTERP_JACOBIAN_FULL)

    # Enable trigger
    if WITH_TRIGGER:
        #ds.eqsys.f_hot.enableIsotropicTrigger(Trigger.TYPE_COLD_ELECTRON_RISE)
        ds.eqsys.f_hot.enableIsotropicTrigger(Trigger.TYPE_TIME, trigger_time=TRIGGER_TIME)
    else:
        ds.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
        ds.eqsys.T_cold.setInitialProfile(1)

    # Solver settings
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_MKL)
    ds.solver.setMaxIterations(100)
    #ds.solver.setVerbose(True)

    ds.solver.tolerance.set(reltol=1e-5)
    ds.solver.setDebug(savejacobian=True, saveresidual=True, timestep=1, iteration=1)

    ds.other.include('fluid', 'scalar')

    ds.setIgnore(['n_i', 'T_cold', 'W_cold', 'n_hot', 'n_cold'])

    # Time step settings
    ds.timestep.setTmax(1e-5)
    ds.timestep.setNt(10000)
    #ds.timestep.setIonization(dt0=1e-8, dtmax=1e-7, tmax=1e-3)

    return ds


def main():
    ds = setup()
    ds.save('settings.h5')

    try:
        if WITH_TRIGGER:
            runiface(ds, 'output_trigger.h5')
        else:
            runiface(ds, 'output_notrigger.h5')
    except: pass

    if WITH_TRIGGER:
        os.rename('petsc_jac', 'petsc_jac_trigger')
        os.rename('residual.mat', 'residual_trigger.mat')
    else:
        os.rename('petsc_jac', 'petsc_jac_notrigger')
        os.rename('residual.mat', 'residual_notrigger.mat')

    return 0


if __name__ == '__main__':
    sys.exit(main())


