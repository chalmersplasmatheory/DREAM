import abc
import os
import pathlib
import sys
from dataclasses import dataclass

import h5py
import numpy as np

from DREAM.DREAMException import DREAMException


GEOMETRY_PACKAGE_MARKER = "dream.stellarator_geometry_package"
GEOMETRY_PACKAGE_V1 = 1
GEOMETRY_PACKAGE_V2 = 2

GROUP_METADATA = "metadata"
GROUP_GRID = "grid"
GROUP_PROFILES = "profiles"
GROUP_SAMPLED = "sampled"
GROUP_BOOZER = "boozer"

REQUIRED_GRID_KEYS = ("rho", "theta", "phi")
REQUIRED_PROFILE_KEYS = ("s", "G", "I", "iota", "psi_T", "B_min", "B_max", "f_passing")
REQUIRED_SAMPLED_KEYS = ("R", "Z", "B", "BdotGradPhi", "Jacobian", "g_tt", "g_tp", "lambda_t", "lambda_p")


def _as_path(filename):
    if filename is None:
        return None

    p = pathlib.Path(filename)
    if p.suffix == "":
        p = p.with_suffix(".h5")

    return p


def _ensure_string(value):
    if value is None:
        return ""
    return str(value)


def _read_string(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _to_numpy(value, dtype=np.float64):
    if value is None:
        return None
    arr = np.asarray(value)
    if dtype is not None and arr.dtype.kind not in ("S", "U", "O"):
        arr = arr.astype(dtype)
    return arr


def _dataset_from_group(group, name, default=None, dtype=np.float64):
    if name not in group:
        return default
    data = group[name][()]
    if isinstance(data, bytes):
        return data.decode("utf-8")
    if np.isscalar(data):
        return data
    return np.asarray(data, dtype=dtype)


def _write_scalar_or_array(group, name, value):
    if isinstance(value, dict):
        sub = group.create_group(name)
        for key, entry in value.items():
            _write_scalar_or_array(sub, key, entry)
        return

    if isinstance(value, (str, pathlib.Path)):
        group.create_dataset(name, data=np.bytes_(str(value)))
        return

    if isinstance(value, (bool, np.bool_)):
        group.create_dataset(name, data=np.uint8(1 if value else 0))
        return

    arr = np.asarray(value)
    if arr.dtype.kind in ("U", "O"):
        group.create_dataset(name, data=np.bytes_(str(value)))
    else:
        group.create_dataset(name, data=arr)


def _read_dict(group):
    out = {}
    for key, value in group.items():
        if isinstance(value, h5py.Group):
            out[key] = _read_dict(value)
        else:
            data = value[()]
            if isinstance(data, bytes):
                out[key] = data.decode("utf-8")
            elif isinstance(data, np.ndarray):
                out[key] = np.asarray(data)
            else:
                out[key] = data
    return out


@dataclass
class StellaratorGeometryPackage:
    metadata: dict
    grid: dict
    profiles: dict
    sampled: dict
    boozer: dict | None = None

    def __post_init__(self):
        self.metadata = dict(self.metadata or {})
        self.grid = {k: _to_numpy(v) for k, v in (self.grid or {}).items()}
        self.profiles = {k: _to_numpy(v) for k, v in (self.profiles or {}).items()}
        self.sampled = {k: _to_numpy(v) for k, v in (self.sampled or {}).items()}
        self.boozer = dict(self.boozer) if self.boozer else None
        self.metadata.setdefault("schema_version", GEOMETRY_PACKAGE_V2 if self.boozer else GEOMETRY_PACKAGE_V1)
        self.validate()

    @property
    def schema_version(self):
        return int(self.metadata.get("schema_version", GEOMETRY_PACKAGE_V1))

    @property
    def nrho(self):
        return int(self.grid["rho"].size)

    @property
    def ntheta(self):
        return int(self.grid["theta"].size)

    @property
    def nphi(self):
        return int(self.grid["phi"].size)

    @property
    def sample_size(self):
        return self.nphi * self.nrho * self.ntheta

    @property
    def a(self):
        return float(self.metadata["minor_radius"])

    @property
    def R0(self):
        return float(self.metadata["major_radius"])

    @property
    def nfp(self):
        return int(self.metadata["nfp"])

    def validate(self):
        for key in REQUIRED_GRID_KEYS:
            if key not in self.grid or self.grid[key] is None:
                raise DREAMException(f"StellaratorGeometryPackage: Missing required grid array '{key}'.")
            if not isinstance(self.grid[key], np.ndarray) or self.grid[key].ndim != 1:
                raise DREAMException(f"StellaratorGeometryPackage: Grid array '{key}' must be 1D.")

        for key in REQUIRED_PROFILE_KEYS:
            if key not in self.profiles or self.profiles[key] is None:
                raise DREAMException(f"StellaratorGeometryPackage: Missing required profile '{key}'.")
            if self.profiles[key].shape != (self.nrho,):
                raise DREAMException(
                    f"StellaratorGeometryPackage: Profile '{key}' must have shape {(self.nrho,)}, got {self.profiles[key].shape}."
                )

        for key in REQUIRED_SAMPLED_KEYS:
            if key not in self.sampled or self.sampled[key] is None:
                raise DREAMException(f"StellaratorGeometryPackage: Missing required sampled array '{key}'.")
            if self.sampled[key].size != self.sample_size:
                raise DREAMException(
                    f"StellaratorGeometryPackage: Sampled array '{key}' must contain {self.sample_size} values, got {self.sampled[key].size}."
                )

        for key in ("major_radius", "minor_radius", "nfp", "source_kind", "provider"):
            if key not in self.metadata:
                raise DREAMException(f"StellaratorGeometryPackage: Missing required metadata field '{key}'.")

    def reshape_sampled(self, key):
        if key not in self.sampled:
            raise DREAMException(f"StellaratorGeometryPackage: No sampled array named '{key}'.")
        return self.sampled[key].reshape((self.nphi, self.nrho, self.ntheta))

    def to_kernel_data(self):
        out = {
            "rho": self.grid["rho"],
            "theta": self.grid["theta"],
            "phi": self.grid["phi"],
            "nfp": self.nfp,
            "R0": self.R0,
            "a": self.a,
        }
        out.update(self.profiles)
        out.update(self.sampled)
        return out

    def write(self, filename):
        path = _as_path(filename)
        if path is None:
            raise DREAMException("StellaratorGeometryPackage: No filename specified when writing package.")

        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "w") as hf:
            hf.attrs[GEOMETRY_PACKAGE_MARKER] = np.uint8(1)
            hf.attrs["schema_version"] = np.uint64(self.schema_version)

            md = hf.create_group(GROUP_METADATA)
            for key, value in self.metadata.items():
                _write_scalar_or_array(md, key, value)

            gg = hf.create_group(GROUP_GRID)
            for key, value in self.grid.items():
                _write_scalar_or_array(gg, key, value)

            gp = hf.create_group(GROUP_PROFILES)
            for key, value in self.profiles.items():
                _write_scalar_or_array(gp, key, value)

            gs = hf.create_group(GROUP_SAMPLED)
            for key, value in self.sampled.items():
                _write_scalar_or_array(gs, key, value)

            if self.boozer:
                gb = hf.create_group(GROUP_BOOZER)
                for key, value in self.boozer.items():
                    _write_scalar_or_array(gb, key, value)

        return path

    @classmethod
    def read(cls, filename):
        path = _as_path(filename)
        if path is None or not path.is_file():
            raise DREAMException(f"StellaratorGeometryPackage: Geometry package '{filename}' does not exist.")

        with h5py.File(path, "r") as hf:
            if GEOMETRY_PACKAGE_MARKER not in hf.attrs:
                raise DREAMException(f"StellaratorGeometryPackage: File '{path}' is not a DREAM geometry package.")

            metadata = _read_dict(hf[GROUP_METADATA])
            metadata["schema_version"] = int(hf.attrs.get("schema_version", GEOMETRY_PACKAGE_V1))
            grid = _read_dict(hf[GROUP_GRID])
            profiles = _read_dict(hf[GROUP_PROFILES])
            sampled = _read_dict(hf[GROUP_SAMPLED])
            boozer = _read_dict(hf[GROUP_BOOZER]) if GROUP_BOOZER in hf else None

        return cls(metadata=metadata, grid=grid, profiles=profiles, sampled=sampled, boozer=boozer)

    @classmethod
    def is_geometry_package(cls, filename):
        path = _as_path(filename)
        if path is None or not path.is_file():
            return False

        try:
            with h5py.File(path, "r") as hf:
                return GEOMETRY_PACKAGE_MARKER in hf.attrs
        except OSError:
            return False

    @classmethod
    def from_legacy_cache(cls, filename, *, source_kind="legacy_cache", provider="legacy_cache", source_path=None):
        path = _as_path(filename)
        if path is None or not path.is_file():
            raise DREAMException(f"StellaratorGeometryPackage: Legacy cache '{filename}' does not exist.")

        with h5py.File(path, "r") as hf:
            rho = np.asarray(hf["rho"][:], dtype=np.float64)
            theta = np.asarray(hf["theta"][:], dtype=np.float64)
            phi = np.asarray(hf["phi"][:], dtype=np.float64)
            a = float(hf["a"][()])
            R0 = float(hf["R0"][()])
            nfp = int(hf["nfp"][()])

            sampled = {
                key: np.asarray(hf[key][:], dtype=np.float64)
                for key in REQUIRED_SAMPLED_KEYS
            }
            profiles = {
                key: np.asarray(hf[key][:], dtype=np.float64)
                for key in ("f_passing", "B_min", "B_max", "G", "I", "iota", "psi_T")
            }

        profiles["s"] = np.square(rho / a) if a > 0 else np.zeros_like(rho)
        metadata = {
            "schema_version": GEOMETRY_PACKAGE_V1,
            "provider": provider,
            "source_kind": source_kind,
            "source_path": _ensure_string(source_path or path),
            "major_radius": R0,
            "minor_radius": a,
            "nfp": nfp,
            "dependency_versions": {},
        }
        grid = {"rho": rho, "theta": theta, "phi": phi}

        return cls(metadata=metadata, grid=grid, profiles=profiles, sampled=sampled, boozer=None)


class GeometryProvider(abc.ABC):
    name = "base"

    def __init__(self, source):
        self.source = str(source)

    @abc.abstractmethod
    def build_package(self, nr, ntheta, nphi, with_boozer=False):
        pass


def import_optional(module_name, package_root_env=None):
    if package_root_env:
        root = os.environ.get(package_root_env)
        if root:
            root = pathlib.Path(root)
            if root.is_dir():
                root_str = str(root)
                if root_str not in sys.path:
                    sys.path.insert(0, root_str)

    try:
        return __import__(module_name, fromlist=["*"])
    except ImportError:
        if package_root_env:
            root = os.environ.get(package_root_env)
            if root:
                root = pathlib.Path(root)
                if root.is_dir():
                    sys.path.insert(0, str(root))
                    try:
                        return __import__(module_name, fromlist=["*"])
                    except ImportError:
                        pass
        raise
