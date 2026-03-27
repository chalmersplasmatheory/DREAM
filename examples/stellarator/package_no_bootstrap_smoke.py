#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

from DREAM.DREAMSettings import DREAMSettings
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.BootstrapCurrent as BootstrapCurrent
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.OhmicCurrent as OhmicCurrent
import DREAM.Settings.RadialGrid as RadialGrid


ROOT = Path(__file__).resolve().parent
DEFAULT_PACKAGE = ROOT / "data" / "stellarator_geometry_v1.h5"
DEFAULT_DREAMI = ROOT.parents[1] / "build-brewpetsc" / "iface" / "dreami"


def build_settings(package: Path, settings_path: Path, output_path: Path, nr: int = 6) -> DREAMSettings:
    ds = DREAMSettings()

    ds.radialgrid.setNr(nr)
    ds.radialgrid.setNtheta(17)
    ds.radialgrid.setNphi(17)
    ds.radialgrid.setStellarator(str(package), provider="package")
    ds.radialgrid.setWallRadius(1.05 * ds.radialgrid.a)

    a = ds.radialgrid.a
    r = np.linspace(0, a, nr)
    rho = r / a

    te = 4.0e3 * (1 - 0.65 * rho**2)
    ne = 5.0e19 * (1 - 0.30 * rho**2)
    ti = 2.5e3 * (1 - 0.40 * rho**2)
    johm = 1 - 0.8 * rho**2

    ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(
        ElectricField.BC_TYPE_PRESCRIBED,
        V_loop_wall_R0=0.0,
        R0=ds.radialgrid.R0,
    )

    ds.eqsys.j_ohm.setInitialProfile(johm, radius=r, Ip0=5.0e4)
    ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)
    ds.eqsys.j_bs.setMode(BootstrapCurrent.BOOTSTRAP_MODE_DISABLED)

    ds.eqsys.T_cold.setPrescribedData(te, radius=r)
    ds.eqsys.T_cold.setType(Temperature.TYPE_PRESCRIBED)

    ds.eqsys.n_i.addIon(
        name="D",
        Z=1,
        iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED,
        n=0.95 * ne,
        r=r,
        T=ti,
    )
    ds.eqsys.n_i.addIon(
        name="C",
        Z=6,
        iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED,
        n=0.05 * ne / 6.0,
        r=r,
        T=ti,
    )

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setMaxIterations(100)
    ds.solver.tolerance.set(reltol=1e-5)

    ds.timestep.setTmax(1.0e-4)
    ds.timestep.setNt(1)

    ds.output.setFilename(str(output_path))
    ds.save(str(settings_path))
    return ds


def run_kernel(dreami: Path, settings_path: Path) -> float:
    t0 = time.perf_counter()
    subprocess.run([str(dreami), str(settings_path)], check=True)
    return time.perf_counter() - t0


def main():
    parser = argparse.ArgumentParser(description="Run a package-backed no-bootstrap DREAM stellarator smoke case.")
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--dreami", type=Path, default=DEFAULT_DREAMI)
    parser.add_argument("--settings", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--nr", type=int, default=6)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    if args.settings is None or args.output is None:
        tmpdir = Path(tempfile.mkdtemp(prefix="dream-stellarator-smoke-"))
        settings_path = tmpdir / "settings.h5"
        output_path = tmpdir / "output.h5"
    else:
        settings_path = args.settings
        output_path = args.output

    build_settings(args.package, settings_path, output_path, nr=args.nr)
    print(f"Settings written to {settings_path}")

    if args.run:
        runtime = run_kernel(args.dreami, settings_path)
        print(f"Kernel execution completed in {runtime:.3f} s")
        print(f"Output written to {output_path}")


if __name__ == "__main__":
    main()
