#!/usr/bin/env python3

import argparse
from pathlib import Path

from DREAM import DREAMSettings
import DREAM.Settings.RadialGrid as RadialGrid


ROOT = Path(__file__).resolve().parent
DEFAULT_SOURCE = ROOT / "data" / "wout_LandremanPaul2021_QA_lowres_reference.nc"
DEFAULT_CACHE = ROOT / "data" / "stellarator_geometry_desc.h5"


def main():
    parser = argparse.ArgumentParser(description="Build a stellarator geometry package through DREAM's DESC-backed frontend.")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--nr", type=int, default=4)
    parser.add_argument("--ntheta-equil", type=int, default=17)
    parser.add_argument("--nphi-equil", type=int, default=17)
    args = parser.parse_args()

    try:
        ds = DREAMSettings()
        ds.radialgrid.setType(RadialGrid.TYPE_STELLARATOR)
        ds.radialgrid.setNr(args.nr - 1)
        ds.radialgrid.setMinorRadius(0.18)
        ds.radialgrid.setWallRadius(0.18)
        ds.radialgrid.setStellarator(
            str(args.source),
            provider="desc",
            cache_filename=str(args.cache),
            write_cache=True,
            nr_equil=args.nr,
            ntheta_equil=args.ntheta_equil,
            nphi_equil=args.nphi_equil,
        )

        print(f"Provider: {ds.radialgrid.stellarator_provider}")
        print(f"Geometry cache: {args.cache}")
        print(f"Major radius: {ds.radialgrid.R0:.8f} m")
        print(f"Minor radius: {ds.radialgrid.a:.8f} m")
        print(f"Sample grid: nrho={ds.radialgrid.rho.size}, ntheta={ds.radialgrid.theta.size}, nphi={ds.radialgrid.phi.size}")
    except Exception as ex:
        print("DESC frontend smoke failed.")
        print(ex)
        print("If DESC is unavailable in this environment, use 'geometry_from_package.py' with a precomputed package instead.")


if __name__ == "__main__":
    main()
