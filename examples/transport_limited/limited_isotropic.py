#!/usr/bin/env python3
#
# Minimal isotropic toroidal test case for the eigenmode-limited
# Rechester-Rosenbluth transport model.
#

from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
sys.path.append(str(REPO_ROOT / "py"))

from DREAM import DREAMOutput, DREAMSettings, runiface
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport


MU0 = 4e-7 * np.pi


def configure_radial_grid(ds, nr, a=0.5, R0=3.0, B0=5.0, Ip_ref=8e5):
    """
    Configure a simple analytic toroidal geometry. The trapped-passing
    boundary is non-trivial in this geometry, which is essential for
    the detrapping limiter to be visible.
    """
    rref = np.linspace(0.0, a, 40)
    psi = -MU0 * Ip_ref * (1.0 - (rref / a) ** 2) * a
    Delta = np.zeros_like(rref)

    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
    ds.radialgrid.setShaping(
        psi=psi,
        rpsi=rref,
        GOverR0=B0,
        kappa=1.0,
        delta=0.0,
        Delta=Delta,
        rDelta=rref
    )
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a)
    ds.radialgrid.setMajorRadius(R0)
    ds.radialgrid.setNr(nr)


def make_settings(
    model,
    dBB=3e-3,
    nr=6,
    np1=80,
    pmax=1.5,
    tmax=1e-6,
    nt=1
):
    """
    Build a minimal isotropic kinetic test case.
    """
    ds = DREAMSettings()

    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL

    configure_radial_grid(ds, nr=nr)

    r = np.linspace(0.0, ds.radialgrid.a, nr)
    t0 = np.array([0.0], dtype=float)
    n_cold = 5e19 * np.ones((1, nr), dtype=float)
    T_cold = 1e3 * np.ones((1, nr), dtype=float)
    E_field = 1e-2 * np.ones((1, nr), dtype=float)

    ds.eqsys.E_field.setPrescribedData(E_field, radius=r, times=t0)
    ds.eqsys.T_cold.setPrescribedData(T_cold, radius=r, times=t0)
    ds.eqsys.n_cold.setPrescribedData(n_cold, radius=r, times=t0)
    ds.eqsys.n_i.addIon(
        name="D",
        Z=1,
        iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED,
        n=5e19
    )
    ds.eqsys.n_re.setDreicer(False)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
    ds.eqsys.n_re.setCompton(False)
    ds.eqsys.n_re.setTritium(False)
    ds.eqsys.n_re.setHottail(False)

    ds.hottailgrid.setEnabled(True)
    ds.runawaygrid.setEnabled(False)
    ds.hottailgrid.setNxi(1)
    ds.hottailgrid.setNp(np1)
    ds.hottailgrid.setPmax(pmax)

    ds.eqsys.f_hot.setInitialProfiles(n0=0.99 * n_cold.flatten(), T0=T_cold.flatten(), rn0=r, rT0=r)
    ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0)
    ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_EXPLICIT)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(
        ad_int=FHot.AD_INTERP_TCDF,
        ad_jac=FHot.AD_INTERP_JACOBIAN_FULL
    )
    ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)

    t_dBB = np.array([0.0, tmax], dtype=float)
    r_dBB = np.linspace(0.0, ds.radialgrid.a, nr + 1)
    dBB_arr = float(dBB) * np.ones((2, nr + 1), dtype=float)

    if model == "rr":
        ds.eqsys.f_hot.transport.setMagneticPerturbation(dBB_arr, t=t_dBB, r=r_dBB)
    elif model == "detrapping":
        ds.eqsys.f_hot.transport.setDetrappingRechesterRosenbluth(
            dBB_arr, t=t_dBB, r=r_dBB
        )
    else:
        raise ValueError(f"Unrecognized model: {model}")

    ds.solver.setType(Solver.LINEAR_IMPLICIT)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.solver.preconditioner.setEnabled(False)

    ds.timestep.setTmax(tmax)
    ds.timestep.setNt(nt)

    ds.other.include("fluid", "hottail/Drr", "hottail/nu_D_f1")

    return ds


def extract_last_drr(do):
    """
    Return the last saved isotropic Drr array with shape (nr_f, np).
    """
    Drr = np.asarray(do.other.hottail.Drr[:], dtype=float)
    if Drr.ndim == 4:
        if Drr.shape[2] != 1:
            raise RuntimeError(f"Expected isotropic Drr with singleton pitch dimension, got {Drr.shape}.")
        return Drr[-1, :, 0, :]
    if Drr.ndim == 3:
        return Drr[-1, :, :]

    raise RuntimeError(f"Unexpected Drr shape: {Drr.shape}.")


def summarize_pair(rr_path, det_path):
    """
    Print a tiny summary showing that the detrapping model reduces Drr.
    """
    do_rr = DREAMOutput(str(rr_path))
    do_det = DREAMOutput(str(det_path))

    D_rr = extract_last_drr(do_rr)
    D_det = extract_last_drr(do_det)
    ratio = np.divide(
        D_det,
        D_rr,
        out=np.full_like(D_det, np.nan),
        where=np.abs(D_rr) > 0
    )

    ir = 1 if D_rr.shape[0] > 1 else 0
    ip = D_rr.shape[1] // 2
    p = float(np.asarray(do_rr.grid.hottail.p[:], dtype=float)[ip])
    r_f = float(np.asarray(do_rr.grid.r_f[:], dtype=float)[ir])

    print("Comparison summary")
    print(f"  radial flux index used : {ir}")
    print(f"  r_f                    : {r_f:.6f} m")
    print(f"  momentum sample        : p = {p:.6f}")
    print(f"  Drr(rr)                : {D_rr[ir, ip]:.6e}")
    print(f"  Drr(detrapping)        : {D_det[ir, ip]:.6e}")
    print(f"  ratio                  : {ratio[ir, ip]:.6e}")
    print(f"  global ratio range     : [{np.nanmin(ratio):.6e}, {np.nanmax(ratio):.6e}]")
    do_rr.close()
    do_det.close()


def run_case(model, args):
    """
    Generate settings, save them, and optionally run DREAM.
    """
    ds = make_settings(
        model=model,
        dBB=args.dBB,
        nr=args.nr,
        np1=args.np,
        pmax=args.pmax,
        tmax=args.tmax,
        nt=args.nt
    )

    tag = "eigenmode" if model == "detrapping" else model
    settings_path = SCRIPT_DIR / f"limited_isotropic_{tag}_settings.h5"
    output_path = SCRIPT_DIR / f"limited_isotropic_{tag}_output.h5"

    ds.output.setFilename(str(output_path))
    ds.save(str(settings_path))

    print(f"[{model}] settings saved to {settings_path}")

    if not args.no_run:
        runiface(ds, str(output_path), quiet=args.quiet)
        print(f"[{model}] output written to {output_path}")

    return settings_path, output_path


def parse_args():
    parser = ArgumentParser(
        description=(
            "Minimal isotropic toroidal RR transport test. By default, runs both "
            "ordinary RR and the trapped-pitch eigenmode model and saves 'hottail/Drr'."
        )
    )
    parser.add_argument(
        "--model",
        choices=["rr", "detrapping", "both"],
        default="both",
        help="Which transport model to run."
    )
    parser.add_argument("--dBB", type=float, default=3e-3, help="Constant magnetic perturbation level.")
    parser.add_argument("--nr", type=int, default=6, help="Number of radial cells.")
    parser.add_argument("--np", type=int, default=80, help="Number of momentum points.")
    parser.add_argument("--pmax", type=float, default=1.5, help="Maximum hot-tail momentum.")
    parser.add_argument("--tmax", type=float, default=1e-6, help="Total simulation time [s].")
    parser.add_argument("--nt", type=int, default=1, help="Number of time steps.")
    parser.add_argument("--no-run", action="store_true", help="Only write settings files, do not run DREAM.")
    parser.add_argument("--quiet", action="store_true", help="Run DREAM quietly.")

    return parser.parse_args()


def main():
    args = parse_args()

    if args.model == "both":
        _, rr_output = run_case("rr", args)
        _, det_output = run_case("detrapping", args)
        if not args.no_run:
            summarize_pair(rr_output, det_output)
    else:
        run_case(args.model, args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
