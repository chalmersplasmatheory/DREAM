#!/usr/bin/env python3

import argparse
from pathlib import Path

from DREAM import DREAMSettings
import DREAM.Settings.RadialGrid as RadialGrid


ROOT = Path(__file__).resolve().parent
DEFAULT_PACKAGE = ROOT / "data" / "stellarator_geometry_v2.h5"


def main():
    parser = argparse.ArgumentParser(description="Load a precomputed DREAM stellarator geometry package and evaluate a sample flux tube.")
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--nr", type=int, default=4)
    args = parser.parse_args()

    ds = DREAMSettings()
    ds.radialgrid.setType(RadialGrid.TYPE_STELLARATOR)
    ds.radialgrid.setNr(args.nr - 1)
    ds.radialgrid.setMinorRadius(0.18)
    ds.radialgrid.setWallRadius(0.18)
    ds.radialgrid.setStellarator(str(args.package), provider="package")

    trace = ds.radialgrid.getStellaratorFluxTubeEvaluator().evaluate(s=0.5, alpha=0.0, nturns=1, npoints=64)

    print(f"Package: {args.package}")
    print(f"Sample grid: nrho={ds.radialgrid.rho.size}, ntheta={ds.radialgrid.theta.size}, nphi={ds.radialgrid.phi.size}")
    print(f"Flux-tube sample points: {trace['R'].size}")
    print(f"Mid-radius iota: {trace['iota']:.8f}")


if __name__ == "__main__":
    main()
