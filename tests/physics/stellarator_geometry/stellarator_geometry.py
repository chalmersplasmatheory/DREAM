#!/usr/bin/env python3

import os
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

import dreamtests

import DREAM.Settings.RadialGrid as RadialGrid
from DREAM.Settings.StellaratorFluxTube import FluxTubeEvaluator
from DREAM.Settings.StellaratorGeometry import StellaratorGeometryPackage
from DREAM.Settings.StellaratorGeometryProviders import DescProvider, VmecJaxProvider


ROOT = pathlib.Path(__file__).resolve().parents[3] / "examples" / "stellarator" / "data"
PACKAGE_V1 = ROOT / "stellarator_geometry_v1.h5"
PACKAGE_V2 = ROOT / "stellarator_geometry_v2.h5"
LEGACY_CACHE = ROOT / "legacy_numeric_stellarator_cache.h5"
WOUT = ROOT / "wout_LandremanPaul2021_QA_lowres_reference.nc"
KERNEL_SMOKE = pathlib.Path(__file__).resolve().parents[3] / "examples" / "stellarator" / "package_no_bootstrap_smoke.py"
DREAMI = pathlib.Path(__file__).resolve().parents[3] / "build-brewpetsc" / "iface" / "dreami"


def _desc_available():
    try:
        DescProvider(str(WOUT)).build_package(nr=4, ntheta=17, nphi=17, with_boozer=False)
        return True
    except Exception:
        return False


def _vmec_jax_available():
    try:
        if "VMEC_JAX_ROOT" not in os.environ:
            return False
        VmecJaxProvider(str(WOUT)).build_package(nr=4, ntheta=17, nphi=17, with_boozer=False)
        return True
    except Exception:
        return False


def _assert(condition, message):
    if not condition:
        raise AssertionError(message)


def test_package_smoke():
    rg = RadialGrid.RadialGrid(ttype=RadialGrid.TYPE_STELLARATOR)
    rg.setNr(3)
    rg.setMinorRadius(0.18)
    rg.setWallRadius(0.18)
    rg.setStellarator(str(PACKAGE_V1), provider="package")

    _assert(rg.B.size == rg.phi.size * rg.rho.size * rg.theta.size, "Unexpected sampled array size.")
    _assert(rg.f_passing.shape == (rg.rho.size,), "Unexpected passing-fraction profile shape.")
    rg.verifySettings()


def test_settings_roundtrip():
    rg = RadialGrid.RadialGrid(ttype=RadialGrid.TYPE_STELLARATOR)
    rg.setNr(3)
    rg.setMinorRadius(0.18)
    rg.setWallRadius(0.18)
    rg.setStellarator(str(PACKAGE_V1), provider="package")

    data = rg.todict(verify=False)
    rg2 = RadialGrid.RadialGrid(ttype=RadialGrid.TYPE_STELLARATOR)
    rg2.fromdict(data)
    rg2.verifySettings()

    _assert(np.allclose(rg.B_min, rg2.B_min), "B_min profile changed after settings roundtrip.")
    _assert(np.allclose(rg.B, rg2.B), "3D B array changed after settings roundtrip.")


def test_cache_readback():
    rg = RadialGrid.RadialGrid(ttype=RadialGrid.TYPE_STELLARATOR)
    rg.setNr(3)
    rg.setMinorRadius(0.18)
    rg.setWallRadius(0.18)
    rg.setStellarator(str(WOUT), provider="package", cache_filename=str(LEGACY_CACHE))
    rg.verifySettings()
    _assert(np.all(np.isfinite(rg.R)), "Legacy cache load produced non-finite R data.")


def test_package_roundtrip():
    pkg = StellaratorGeometryPackage.read(PACKAGE_V1)
    with tempfile.TemporaryDirectory() as td:
        out = pathlib.Path(td) / "roundtrip_package.h5"
        pkg.write(out)
        reloaded = StellaratorGeometryPackage.read(out)

    for key in pkg.sampled:
        _assert(np.allclose(pkg.sampled[key], reloaded.sampled[key]), f"Sampled array '{key}' changed after package roundtrip.")


def test_flux_tube_evaluator():
    ft = FluxTubeEvaluator(PACKAGE_V2)
    out = ft.evaluate(s=0.5, alpha=0.1, nturns=1, npoints=64)
    _assert(out["R"].shape == (64,), "Flux-tube R trace has the wrong shape.")
    _assert(np.all(np.isfinite(out["|B|"])), "Flux-tube trace contains non-finite B values.")


def test_kernel_smoke():
    if not DREAMI.is_file():
        return "skip"

    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        settings = td / "settings.h5"
        output = td / "output.h5"
        subprocess.run(
            [
                sys.executable,
                str(KERNEL_SMOKE),
                "--package",
                str(PACKAGE_V1),
                "--dreami",
                str(DREAMI),
                "--settings",
                str(settings),
                "--output",
                str(output),
                "--run",
            ],
            check=True,
            env=dict(os.environ, PYTHONPATH=str(pathlib.Path(__file__).resolve().parents[3] / "py")),
        )
        _assert(output.is_file(), "Kernel smoke test did not produce an output file.")


def test_desc_optional():
    if not _desc_available():
        return "skip"

    pkg = DescProvider(str(WOUT)).build_package(nr=4, ntheta=17, nphi=17, with_boozer=False)
    ref = StellaratorGeometryPackage.read(PACKAGE_V1)
    _assert(np.allclose(pkg.profiles["iota"], ref.profiles["iota"], rtol=1e-4, atol=1e-6), "DESC iota profile drifted from the reference package.")


def test_vmec_optional():
    if not _vmec_jax_available():
        return "skip"

    pkg = VmecJaxProvider(str(WOUT)).build_package(nr=4, ntheta=17, nphi=17, with_boozer=False)
    ref = StellaratorGeometryPackage.read(PACKAGE_V1)
    iota_on_ref = np.interp(ref.profiles["s"], pkg.profiles["s"], pkg.profiles["iota"])
    _assert(np.allclose(iota_on_ref, ref.profiles["iota"], rtol=3e-4, atol=1e-6), "vmec_jax iota profile drifted from the reference package.")
    _assert(pkg.boozer is not None, "vmec_jax provider did not attach a spectral block.")


def run(args):
    tests = [
        ("package_smoke", test_package_smoke),
        ("settings_roundtrip", test_settings_roundtrip),
        ("cache_readback", test_cache_readback),
        ("package_roundtrip", test_package_roundtrip),
        ("flux_tube_evaluator", test_flux_tube_evaluator),
        ("kernel_smoke", test_kernel_smoke),
        ("desc_optional", test_desc_optional),
        ("vmec_optional", test_vmec_optional),
    ]

    success = True
    skipped = []

    for name, test in tests:
        try:
            result = test()
            if result == "skip":
                skipped.append(name)
                continue
            dreamtests.print_ok(name)
        except Exception as ex:
            dreamtests.print_error(f"{name}: {ex}")
            success = False

    if skipped and args.get("verbose", False):
        print("Skipped optional tests: {}".format(", ".join(skipped)))

    return success
