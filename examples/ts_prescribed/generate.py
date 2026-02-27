#!/usr/bin/env python3
"""
Example: Prescribed time stepper (TYPE_PRESCRIBED) to resolve initial transients
in a self-consistent current quench scenario.

This script writes two settings files:
  - init_settings.h5      (initial transient resolution)
  - restart_settings.h5   (restart with self-consistent E and optional T_cold)
"""

import numpy as np
import sys
from pathlib import Path

try:
    from DREAM.DREAMSettings import DREAMSettings
except ModuleNotFoundError:
    sys.path.append(str(Path(__file__).resolve().parents[2] / "py"))
    from DREAM.DREAMSettings import DREAMSettings

import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.ElectricField as Efield
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.ColdElectronTemperature as T_cold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TimeStepper as TimeStepper

INIT_SETTINGS_NAME = "init_settings.h5"
INIT_OUTPUT_NAME = "output_init.h5"
RESTART_SETTINGS_NAME = "restart_settings.h5"
RESTART_OUTPUT_NAME = "output.h5"


def quadratic_time_grid(tmax: float, n: int) -> np.ndarray:
    """
    Time grid clustered near t=0, ending exactly at tmax.
    """
    x = np.linspace(0.0, 1.0, n)
    return (x * x) * tmax


def piecewise_linear_time_grid(tmax: float) -> np.ndarray:
    """
    Resolve an initial transient (0..1e-5) with fine steps, then slower CQ.
    """
    t1 = np.linspace(0.0, 1e-5, 5)
    t2 = np.linspace(2e-5, tmax, 20)
    t = np.concatenate((t1, t2))
    return t


def configure_collisions(ds: DREAMSettings) -> None:
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT


def configure_solver(ds: DREAMSettings) -> None:
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.tolerance.set(reltol=1e-4)
    ds.solver.setMaxIterations(maxiter=100)
    ds.solver.setVerbose(True)


def configure_hottail(ds: DREAMSettings, *, enabled: bool, nxi: int, np: int, pmax: float, T0: float) -> None:
    if not enabled:
        ds.hottailgrid.setEnabled(False)
        ds.runawaygrid.setEnabled(False)
        return

    ds.hottailgrid.setNxi(nxi)
    ds.hottailgrid.setNp(np)
    ds.hottailgrid.setPmax(pmax)

    nfree0, rn0 = ds.eqsys.n_i.getFreeElectronDensity()
    ds.eqsys.f_hot.setInitialProfiles(rn0=rn0, n0=nfree0, rT0=0, T0=T0)
    ds.eqsys.f_hot.setBoundaryCondition(bc=FHot.BC_F_0)
    ds.eqsys.f_hot.enableIonJacobian(False)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int=FHot.AD_INTERP_TCDF)

    # Disable runaway grid for this example
    ds.runawaygrid.setEnabled(False)


# ---------------------------------------------------------------------------
# Build settings
# ---------------------------------------------------------------------------
def build_initial_settings(
    *,
    B0: float,
    E_initial: float,
    T_initial: float,
    Nr: int,
    radius: np.ndarray,
    radius_wall: float,
    tmax_init: float,
    n_init_times: int,
    nxi: int,
    np: int,
    pmax: float,
    hottail_enabled: bool,
) -> DREAMSettings:
    ds = DREAMSettings()

    configure_collisions(ds)

    # Radial grid
    ds.radialgrid.setB0(B0)
    ds.radialgrid.setMinorRadius(float(radius[-1]))
    ds.radialgrid.setWallRadius(radius_wall)
    ds.radialgrid.setNr(Nr)

    # Prescribed time stepper: cluster points near t=0
    ds.timestep.setType(TimeStepper.TYPE_PRESCRIBED)
    ds.timestep.setTimes(quadratic_time_grid(tmax_init, n_init_times))
    ds.timestep.setVerbose(True)

    # Ions
    ds.eqsys.n_i.addIon(name="D",  Z=1,  iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=1e20)
    ds.eqsys.n_i.addIon(name="Ar", Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL,      n=1e20)

    # Prescribed E and T
    times_param = np.array([0.0])
    radius_param = np.array([0.0, float(radius[-1])])

    efield = E_initial * np.ones((times_param.size, radius_param.size))
    ds.eqsys.E_field.setPrescribedData(efield=efield, times=times_param, radius=radius_param)

    temperature = T_initial * np.ones((times_param.size, radius_param.size))
    ds.eqsys.T_cold.setPrescribedData(temperature=temperature, times=times_param, radius=radius_param)

    configure_hottail(ds, enabled=hottail_enabled, nxi=nxi, np=np, pmax=pmax, T0=T_initial)

    configure_solver(ds)

    ds.other.include("fluid")

    # Output
    ds.output.setFilename(INIT_OUTPUT_NAME)
    return ds


def build_restart_settings(
    ds_init: DREAMSettings,
    *,
    E_wall: float,
    tmax_restart: float,
    T_selfconsistent: bool,
) -> DREAMSettings:
    ds = DREAMSettings(ds_init, chain=True)

    # Switch to self-consistent electric field
    ds.eqsys.E_field.setType(Efield.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(
        bctype=Efield.BC_TYPE_PRESCRIBED,
        inverse_wall_time=0,
        V_loop_wall_R0=E_wall * 2*np.pi
    )

    # Optional self-consistent cold temperature evolution
    if T_selfconsistent:
        ds.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)

    # Prescribed time grid for restart (fine early, then coarser)
    ds.timestep.setType(TimeStepper.TYPE_PRESCRIBED)
    ds.timestep.setTimes(piecewise_linear_time_grid(tmax_restart))

    ds.output.setFilename(RESTART_OUTPUT_NAME)
    return ds


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # Simulation toggles
    T_SELFCONSISTENT = True
    HOTTAIL_ENABLED = True

    # Physical / numerical parameters
    B0 = 5.0                 # T
    E_INITIAL = 60.0         # V/m
    E_WALL = 0.0             # V/m (boundary)
    T_INITIAL = 6.0          # eV

    NR = 4
    NP = 40
    NXI = 7
    PMAX = 0.03              # m_e*c

    RADIUS = np.array([0.0, 1.0])   # minor radius span (m)
    R_WALL = 1.5                    # wall radius (m)

    TMAX_INIT = 1e-3
    N_INIT_TIMES = 20

    TMAX_RESTART = 5e-3

    # Build + save initial stage
    ds_init = build_initial_settings(
        B0=B0,
        E_initial=E_INITIAL,
        T_initial=T_INITIAL,
        Nr=NR,
        radius=RADIUS,
        radius_wall=R_WALL,
        tmax_init=TMAX_INIT,
        n_init_times=N_INIT_TIMES,
        nxi=NXI,
        np=NP,
        pmax=PMAX,
        hottail_enabled=HOTTAIL_ENABLED,
    )
    ds_init.save(INIT_SETTINGS_NAME)

    ds_restart = build_restart_settings(
        ds_init,
        E_wall=E_WALL,
        tmax_restart=TMAX_RESTART,
        T_selfconsistent=T_SELFCONSISTENT,
    )
    ds_restart.save(RESTART_SETTINGS_NAME)


if __name__ == "__main__":
    main()
