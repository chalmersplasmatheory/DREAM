#!/usr/bin/env python3

import argparse
from pathlib import Path

from DREAM import DREAMSettings
import DREAM.Settings.RadialGrid as RadialGrid


ROOT = Path(__file__).resolve().parent
DEFAULT_WOUT = ROOT / "data" / "wout_LandremanPaul2021_QA_lowres_reference.nc"
DEFAULT_CACHE = ROOT / "data" / "stellarator_geometry_v2.h5"


def main():
    parser = argparse.ArgumentParser(description="Build a stellarator geometry package from a VMEC wout file through the DREAM frontend.")
    parser.add_argument("--wout", type=Path, default=DEFAULT_WOUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--nr", type=int, default=4)
    parser.add_argument("--ntheta-equil", type=int, default=17)
    parser.add_argument("--nphi-equil", type=int, default=17)
    parser.add_argument("--with-boozer", action="store_true")
    args = parser.parse_args()

    try:
        ds = DREAMSettings()
        ds.radialgrid.setType(RadialGrid.TYPE_STELLARATOR)
        ds.radialgrid.setNr(args.nr - 1)
        ds.radialgrid.setMinorRadius(0.18)
        ds.radialgrid.setWallRadius(0.18)
        ds.radialgrid.setStellarator(
            str(args.wout),
            provider="vmec_jax",
            cache_filename=str(args.cache),
            write_cache=True,
            nr_equil=args.nr,
            ntheta_equil=args.ntheta_equil,
            nphi_equil=args.nphi_equil,
            with_boozer=args.with_boozer,
        )

        ft = ds.radialgrid.getStellaratorFluxTubeEvaluator()
        trace = ft.evaluate(s=0.5, alpha=0.0, nturns=1, npoints=64)

        print(f"Provider: {ds.radialgrid.stellarator_provider}")
        print(f"Geometry cache: {args.cache}")
        print(f"Sample grid: nrho={ds.radialgrid.rho.size}, ntheta={ds.radialgrid.theta.size}, nphi={ds.radialgrid.phi.size}")
        print(f"Flux-tube sample points: {trace['R'].size}")
        print(f"Mid-radius iota: {trace['iota']:.8f}")
    except Exception as ex:
        print("VMEC frontend smoke failed.")
        print(ex)
        print("If DESC is unavailable in this environment, use 'geometry_from_package.py' with a precomputed package instead.")


if __name__ == "__main__":
    main()
