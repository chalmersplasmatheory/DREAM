#!/usr/bin/env python3
"""
h5file_to_imas.py

Map one DREAM HDF5 output file into selected IMAS IDSs with IMAS-Python.


Current mapping overview
------------------------
The detailed node-by-node mapping is written to the Markdown report generated
next to the output. The main mapping entry points are:

| IMAS IDS | Main DREAM sources | Mapping functions |
| --- | --- | --- |
| plasma_profiles | /eqsys/T_cold, /n_tot or /n_cold, /n_i, /j_tot, /j_ohm, /E_field, /other/fluid/Zeff, /conductivity | map_plasma_profiles(), fill_ion_profiles_1d(), fill_1d_grid() |
| runaway_electrons | /eqsys/n_re, /j_re, /other/fluid/runawayRate, /gamma*, /pCrit*, /E* | map_runaway_electrons(), fill_1d_grid() |
| spi | /settings/eqsys/spi/init/*, /settings/eqsys/n_i/SPIMolarFraction, /eqsys/x_p, /v_p, /Y_p | map_spi(), map_spi_species(), spi_injector_groups() |
| equilibrium | /eqsys/I_p, /psi_p, /j_tot, /grid/R0, /grid/r, /grid/eq/*, /grid/geometry/* | map_equilibrium(), fill_equilibrium_profiles_2d(), boundary_outline() |
| radiation | /other/fluid/Tcold_radiation, /grid/t, /grid/VpVol, /grid/dr, /grid/R0 | map_radiation() |
| summary | volume/current reductions of plasma/runaway/equilibrium data, /code/* provenance | map_summary(), fill_code_metadata(), volume_integral_trace(), current_from_j_trace() |

How to add a new mapped field
-----------------------------
1. Find the relevant `map_<ids>()` function above.
2. Read the DREAM dataset with `dream.arr("/path", report)` or `dream.scalar()`.
3. Align time-dependent arrays with `time_aligned()` if needed.
4. Fill the IMAS node with `fill_ids_field(root, "imas/node", value, report, ids_name, "DREAM source or derivation")`.
5. Put the derivation in the last argument. This text is what appears in the
   mapping report, so keep it precise enough for DREAM and IMAS experts to audit.

`fill_ids_field()` is intentionally defensive: if an IMAS node is absent in the
installed Data Dictionary, or the assignment fails, the converter records a
skipped mapping in the report instead of crashing.

Units policy
------------
Values are kept in DREAM units where these match IMAS nodes. DREAM T_cold is eV,
densities are m^-3, current densities are A m^-2, E_field is V m^-1, time is
seconds, and r/R0 are metres. Explicit conversions are noted in report sources
for the affected fields.

Requirements
------------
    pip install h5py numpy imas-python

For writing to an IMAS HDF5 data entry you also need an IMAS-Core installation
with the hdf5 backend available. Without IMAS-Core, use `--backend netcdf`.

Examples
--------
Write to a local IMAS HDF5 data-entry directory:

    python dream_to_imas.py output_CQ_S6.h5 --uri 'imas:hdf5?path=./imas_dream'

Write to a portable netCDF file:

    python dream_to_imas.py output_CQ_S6.h5 --uri dream_imas.nc

Only build IDSs and produce a report, without writing:

    python dream_to_imas.py output_CQ_S6.h5 --dry-run
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import h5py
import numpy as np

try:
    import imas
except Exception as exc:  # pragma: no cover - depends on user environment
    imas = None  # type: ignore[assignment]
    IMAS_IMPORT_ERROR = exc
else:
    IMAS_IMPORT_ERROR = None

# ------------------------ normalization factors---------------------------
psi_cocos = 2.0*np.pi
phi_sign   = -1.0  # positive is anti-clockwise in IMAS, ITER currents are negative in IMAS
c_light    = 299792458.0
m_electron = 9.10938356e-31
p_norm = m_electron * c_light

DEFAULT_IDS = ["plasma_profiles", "runaway_electrons", "spi", "equilibrium", "radiation", "summary"]


# ----------------------------- reporting ---------------------------------

@dataclass
class MappingReport:
    source_file: str
    written_uri: Optional[str] = None
    dd_version_requested: Optional[str] = None
    created_ids: list[str] = field(default_factory=list)
    set_nodes: dict[tuple[str, str, str], int] = field(default_factory=dict)
    skipped_nodes: dict[tuple[str, str, str, str], int] = field(default_factory=dict)
    missing_sources: dict[str, int] = field(default_factory=dict)
    warnings: dict[str, int] = field(default_factory=dict)
    timings: dict[str, float] = field(default_factory=dict)
    object_paths: dict[int, str] = field(default_factory=dict, repr=False)

    def bind(self, obj: Any, node_path: str) -> None:
        self.object_paths[id(obj)] = normalize_node_path(node_path)

    def full_node_path(self, ids: str, target: str, root: Any | None = None) -> str:
        base = self.object_paths.get(id(root), ids) if root is not None else ids
        target = normalize_report_text(target)
        if not target:
            return base
        return normalize_node_path(f"{base}.{target.replace('/', '.')}")

    def ok(self, ids: str, target: str, source: str, units: str = "", root: Any | None = None) -> None:
        node = self.full_node_path(ids, target, root)
        key = (
            node,
            units or imas_units_for_node(node),
            normalize_report_text(source),
        )
        self.set_nodes[key] = self.set_nodes.get(key, 0) + 1

    def skip(
        self,
        ids: str,
        target: str,
        reason: str,
        source: str = "",
        units: str = "",
        root: Any | None = None,
    ) -> None:
        node = self.full_node_path(ids, target, root)
        key = (
            node,
            units or imas_units_for_node(node),
            normalize_report_text(source),
            normalize_report_text(reason),
        )
        self.skipped_nodes[key] = self.skipped_nodes.get(key, 0) + 1

    def missing(self, source: str) -> None:
        key = normalize_report_text(source)
        self.missing_sources[key] = self.missing_sources.get(key, 0) + 1

    def warn(self, message: str) -> None:
        self.warnings[message] = self.warnings.get(message, 0) + 1

    def add_timing(self, key: str, seconds: float) -> None:
        self.timings[key] = self.timings.get(key, 0.0) + seconds

    def write(self, path: Path) -> None:
        if path.suffix.lower() == ".json":
            self.write_json(path)
        else:
            self.write_markdown(path)

    def write_json(self, path: Path) -> None:
        payload = {
            "source_file": self.source_file,
            "written_uri": self.written_uri,
            "dd_version_requested": self.dd_version_requested,
            "created_ids": self.created_ids,
            "set_nodes": [
                {"node": node, "units": units, "source": source, "count": count}
                for (node, units, source), count in sorted_report_items(self.set_nodes)
            ],
            "skipped_nodes": [
                {"node": node, "units": units, "source": source, "reason": reason, "count": count}
                for (node, units, source, reason), count in sorted_report_items(self.skipped_nodes)
            ],
            "missing_sources": [
                {"source": source, "count": count}
                for source, count in sorted(self.missing_sources.items(), key=lambda item: (-item[1], item[0]))
            ],
            "warnings": [
                {"message": message, "count": count}
                for message, count in sorted(self.warnings.items(), key=lambda item: (-item[1], item[0]))
            ],
            "timings": [
                {"section": key, "seconds": seconds}
                for key, seconds in sorted(self.timings.items(), key=lambda item: item[0])
            ],
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def write_markdown(self, path: Path) -> None:
        lines = [
            "# DREAM to IMAS Mapping Report",
            "",
            "## Summary",
            "",
            f"- Source file: `{self.source_file}`",
            f"- Written URI: `{self.written_uri or '(dry run)'}`",
            f"- Data Dictionary version requested: `{self.dd_version_requested or '(environment default)'}`",
            f"- Created IDSs: {', '.join(f'`{ids}`' for ids in self.created_ids) if self.created_ids else '(none)'}",
            f"- Set mappings: {len(self.set_nodes)} unique, {sum(self.set_nodes.values())} total calls",
            f"- Skipped mappings: {len(self.skipped_nodes)} unique, {sum(self.skipped_nodes.values())} total calls",
            f"- Missing sources: {len(self.missing_sources)} unique, {sum(self.missing_sources.values())} total calls",
            f"- Warnings: {len(self.warnings)} unique, {sum(self.warnings.values())} total calls",
            "",
        ]

        lines.extend(markdown_table(
            "Set Mappings",
            ["IMAS Node", "Units", "DREAM Source / Derivation", "Count"],
            [(node, units, source, str(count)) for (node, units, source), count in sorted_report_items(self.set_nodes)],
        ))
        lines.extend(markdown_table(
            "Skipped Mappings",
            ["IMAS Node", "Units", "DREAM Source / Derivation", "Reason", "Count"],
            [
                (node, units, source, reason, str(count))
                for (node, units, source, reason), count in sorted_report_items(self.skipped_nodes)
            ],
        ))
        lines.extend(markdown_table(
            "Missing Sources",
            ["Source", "Count"],
            [(source, str(count)) for source, count in sorted(self.missing_sources.items(), key=lambda item: (-item[1], item[0]))],
        ))
        lines.extend(markdown_table(
            "Warnings",
            ["Message", "Count"],
            [(message, str(count)) for message, count in sorted(self.warnings.items(), key=lambda item: (-item[1], item[0]))],
        ))
        lines.extend(markdown_table(
            "Timings",
            ["Section", "Seconds"],
            [(key, f"{seconds:.3f}") for key, seconds in sorted(self.timings.items(), key=lambda item: item[1], reverse=True)],
        ))

        path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@dataclass
class ConversionResult:
    source_file: str
    ids_list: list[Any]
    report: MappingReport


def normalize_report_text(text: str) -> str:
    text = str(text)
    text = re.sub(r"\[:,\d+,:\]", "[:,*,:]", text)
    text = re.sub(r"\[:,\d+\]", "[:,*]", text)
    text = re.sub(r"\[\d+\]", "[*]", text)
    text = re.sub(r"/\d+(?=/|$)", "/*", text)
    return text


def normalize_node_path(node_path: str) -> str:
    node_path = normalize_report_text(node_path)
    node_path = node_path.replace("/", ".")
    node_path = re.sub(r"\.\*?(?=\.|$)", lambda match: ".*" if match.group(0).startswith(".*") else ".", node_path)
    node_path = re.sub(r"\.+", ".", node_path)
    return node_path.strip(".")


IMAS_UNIT_FALLBACKS = [
    ("*.name", "n/a"),
    ("*.label", "n/a"),
    ("*.description", "n/a"),
    ("*.index", "1"),
    ("*.code.*", "n/a"),
    ("*.time", "s"),
    ("*.global_quantities.ip", "A"),
    ("*.global_quantities.ip.value", "A"),
    ("*.global_quantities.r0", "m"),
    ("*.global_quantities.r0.value", "m"),
    ("*.global_quantities.b0", "T"),
    ("*.global_quantities.b0.value", "T"),
    ("*.global_quantities.current", "A"),
    ("*.global_quantities.current_phi", "A"),
    ("*.global_quantities.current_ohm.value", "A"),
    ("*.global_quantities.energy*.value", "J"),
    ("*.global_quantities.psi*.value", "Wb"),
    ("*.vacuum_toroidal_field.r0", "m"),
    ("*.vacuum_toroidal_field.b0", "T"),
    ("*.profiles_1d.grid.rho_tor_norm", "1"),
    ("*.profiles_1d.grid.rho_tor", "m"),
    ("*.profiles_1d.grid.rho_pol_norm", "1"),
    ("*.profiles_1d.grid.volume", "m^3"),
    ("*.profiles_1d.grid.area", "m^2"),
    ("*.profiles_1d.electrons.temperature", "eV"),
    ("*.profiles_1d.electrons.density*", "m^-3"),
    ("*.profiles_1d.ion*.density*", "m^-3"),
    ("*.profiles_1d.neutral*.density*", "m^-3"),
    ("*.profiles_1d.*.density", "m^-3"),
    ("*.profiles_1d.current_density*", "A m^-2"),
    ("*.profiles_1d.j_*", "A m^-2"),
    ("*.profiles_1d.e_field*", "V m^-1"),
    ("*.profiles_1d.conductivity*", "Ohm^-1 m^-1"),
    ("*.profiles_1d.zeff", "1"),
    ("*.profiles_1d.rho_tor_norm", "1"),
    ("*.profiles_1d.rho_tor", "m"),
    ("*.profiles_1d.phi_tor", "Wb"),
    ("*.profiles_1d.r_outboard", "m"),
    ("*.profiles_1d.b_field*", "T"),
    ("*.profiles_1d.trapped_fraction", "1"),
    ("*.profiles_1d.gm1", "1"),
    ("*.profiles_1d.gm5", "T^2"),
    ("*.profiles_1d.momentum_critical*", "kg m s^-1"),
    ("*.profiles_1d.ddensity_dt*", "m^-3 s^-1"),
    ("*.boundary.outline.r", "m"),
    ("*.boundary.outline.z", "m"),
    ("*.boundary.geometric_axis.r", "m"),
    ("*.boundary.geometric_axis.z", "m"),
    ("*.magnetic_axis.r", "m"),
    ("*.magnetic_axis.z", "m"),
    ("*.profiles_2d.grid.dim1", "Wb"),
    ("*.profiles_2d.grid.dim2", "rad"),
    ("*.profiles_2d.r", "m"),
    ("*.profiles_2d.z", "m"),
    ("*.profiles_2d.psi", "Wb"),
    ("summary.volume_average.t_e.value", "eV"),
    ("summary.volume_average.n_e.value", "m^-3"),
    ("summary.volume_average.zeff.value", "1"),
    ("summary.local.*.t_e.value", "eV"),
    ("summary.local.*.n_e.value", "m^-3"),
    ("summary.local.*.n_re.value", "m^-3"),
    ("summary.local.*.j_ohm.value", "A m^-2"),
    ("summary.local.*.j_re.value", "A m^-2"),
    ("summary.local.*.e_field.value", "V m^-1"),
    ("summary.local.*.zeff.value", "1"),
    ("radiation.process.identifier.index", "1"),
    ("radiation.process.identifier.name", "n/a"),
    ("radiation.process.identifier.description", "n/a"),
    ("radiation.process.global_quantities.time", "s"),
    ("radiation.process.global_quantities.inside_vessel.power*", "W"),
    ("radiation.process.profiles_1d.time", "s"),
    ("radiation.process.profiles_1d.grid.rho_tor_norm", "1"),
    ("radiation.process.profiles_1d.grid.rho_tor", "m"),
    ("radiation.process.profiles_1d.electrons.emissivity", "W m^-3"),
    ("radiation.process.profiles_1d.electrons.power_inside", "W"),
    ("spi.injector.pellet.core.atoms_n", "1"),
    ("spi.injector.pellet.core.species.density", "m^-3"),
    ("spi.injector.pellet.core.species.z_n", "1"),
    ("spi.injector.pellet.core.species.a", "1"),
    ("spi.injector.fragment.position.r", "m"),
    ("spi.injector.fragment.position.phi", "rad"),
    ("spi.injector.fragment.position.z", "m"),
    ("spi.injector.fragment.volume", "m^3"),
    ("spi.injector.fragment.velocity_r", "m s^-1"),
    ("spi.injector.fragment.velocity_z", "m s^-1"),
    ("spi.injector.fragment.velocity_phi", "m s^-1"),
    ("spi.injector.shattering_position.r", "m"),
    ("spi.injector.shattering_position.phi", "rad"),
    ("spi.injector.shattering_position.z", "m"),
    ("spi.injector.shatter_cone.origin.r", "m"),
    ("spi.injector.shatter_cone.origin.phi", "rad"),
    ("spi.injector.shatter_cone.origin.z", "m"),
    ("spi.injector.velocity_mass_centre_fragments_r", "m s^-1"),
    ("spi.injector.velocity_mass_centre_fragments_z", "m s^-1"),
    ("spi.injector.velocity_mass_centre_fragments_phi", "m s^-1"),
    ("*.z_n", "1"),
    ("*.a", "1"),
    ("*.density", "m^-3"),
    ("*.temperature", "eV"),
    ("*.current_density", "A m^-2"),
    ("*.current", "A"),
    ("*.e_field*", "V m^-1"),
    ("*.psi*", "Wb"),
    ("*.r", "m"),
    ("*.z", "m"),
    ("*.volume", "m^3"),
]


def imas_units_for_node(node: str) -> str:
    for pattern, units in IMAS_UNIT_FALLBACKS:
        if fnmatch.fnmatchcase(node, pattern):
            return units
    return ""


def imas_units_from_node(node: Any) -> str:
    for attr in ("units", "unit"):
        try:
            value = getattr(node, attr)
        except Exception:
            value = None
        if value:
            return str(value)

    for attr in ("metadata", "_metadata"):
        try:
            metadata = getattr(node, attr)
        except Exception:
            metadata = None
        units = metadata_units(metadata)
        if units:
            return units

    return ""


def metadata_units(metadata: Any) -> str:
    if metadata is None:
        return ""
    if isinstance(metadata, dict):
        for key in ("units", "unit", "documentation_unit"):
            value = metadata.get(key)
            if value:
                return str(value)
    for attr in ("units", "unit", "documentation_unit"):
        try:
            value = getattr(metadata, attr)
        except Exception:
            value = None
        if value:
            return str(value)
    return ""


def sorted_report_items(items: dict[tuple[str, ...], int]) -> list[tuple[tuple[str, ...], int]]:
    return sorted(items.items(), key=lambda item: (-item[1], item[0]))


def escape_markdown_cell(value: Any) -> str:
    text = str(value)
    text = text.replace("\\", "\\\\").replace("|", "\\|")
    text = text.replace("\n", "<br>")
    return text


def markdown_table(title: str, headers: list[str], rows: list[tuple[Any, ...]], limit: int = 250) -> list[str]:
    lines = [f"## {title}", ""]
    if not rows:
        lines.extend(["No entries.", ""])
        return lines

    visible_rows = rows[:limit]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in visible_rows:
        lines.append("| " + " | ".join(escape_markdown_cell(value) for value in row) + " |")
    if len(rows) > limit:
        lines.extend(["", f"Showing the first {limit} rows of {len(rows)} unique entries.", ""])
    else:
        lines.append("")
    return lines


# ----------------------------- HDF5 helpers -------------------------------

class DreamH5:
    def __init__(self, filename: str | Path):
        self.filename = str(filename)
        self.h5 = h5py.File(self.filename, "r")
        self._cache: dict[str, Any] = {}

    def close(self) -> None:
        self.h5.close()

    def exists(self, path: str) -> bool:
        return path in self.h5

    def arr(self, path: str, report: MappingReport | None = None) -> Optional[np.ndarray]:
        if path in self._cache:
            return self._cache[path]
        if path not in self.h5:
            if report is not None:
                report.missing(path)
            return None
        data = self.h5[path][()]
        value = decode_h5_value(data)
        self._cache[path] = value
        return value

    def scalar(self, path: str, default: float | None = None) -> Optional[float]:
        value = self.arr(path)
        if value is None:
            return default
        arr = np.asarray(value)
        if arr.size == 0:
            return default
        return float(arr.reshape(-1)[0])


def decode_h5_value(value: Any) -> Any:
    """Decode bytes/character arrays that appear in DREAM HDF5 files."""
    arr = np.asarray(value)
    if arr.dtype.kind in {"S", "U"}:
        try:
            if arr.ndim == 0:
                return arr.item().decode() if isinstance(arr.item(), bytes) else str(arr.item())
            # DREAM often stores strings as arrays of single characters.
            chars = []
            for x in arr.reshape(-1):
                if isinstance(x, bytes):
                    chars.append(x.decode(errors="ignore"))
                else:
                    chars.append(str(x))
            return "".join(chars).rstrip("\x00")
        except Exception:
            return arr.astype(str)
    return arr


def flatten_1d(x: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if x is None:
        return None
    a = np.asarray(x)
    if a.size == 0:
        return a.astype(float)
    return np.asarray(a, dtype=float).reshape(-1)


def scale_optional(data: Optional[np.ndarray], factor: float) -> Optional[np.ndarray]:
    if data is None:
        return None
    return np.asarray(data, dtype=float) * factor


def as_text(value: Any) -> Optional[str]:
    """Convert DREAM scalar/string datasets to a compact Python string."""
    if value is None:
        return None
    if isinstance(value, str):
        text = value
    elif isinstance(value, bytes):
        text = value.decode(errors="ignore")
    else:
        arr = np.asarray(value)
        if arr.size == 0:
            return None
        if arr.dtype.kind in {"S", "U", "O"}:
            text = "".join(str(item.decode(errors="ignore") if isinstance(item, bytes) else item) for item in arr.reshape(-1))
        else:
            text = str(arr.reshape(-1)[0])
    text = text.strip().rstrip("\x00")
    return text or None


def dream_text(dream: DreamH5, path: str, report: MappingReport | None = None) -> Optional[str]:
    return as_text(dream.arr(path, report))


def xml_text(value: str) -> str:
    """Escape a string for the small XML fragment stored in code/parameters."""
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def dream_code_parameters_xml(values: dict[str, str]) -> str:
    lines = ["<dream_code_metadata>"]
    for key, value in values.items():
        lines.append(f"  <{key}>{xml_text(value)}</{key}>")
    lines.append("</dream_code_metadata>")
    return "\n".join(lines)


def time_aligned(data: Optional[np.ndarray], nt: int) -> Optional[np.ndarray]:
    """Return data aligned to a time dimension, padding edge values if needed."""
    if data is None:
        return None
    a = np.asarray(data, dtype=float)
    if a.ndim == 0:
        return np.full((nt,), float(a))
    if a.shape[0] == nt:
        return a
    if a.shape[0] == nt - 1:
        pad = np.take(a, [-1], axis=0)
        return np.concatenate([a, pad], axis=0)
    if a.shape[0] > nt:
        return a[:nt]
    if a.shape[0] < nt and a.shape[0] > 0:
        reps = [1] * a.ndim
        reps[0] = nt - a.shape[0]
        pad = np.repeat(np.take(a, [-1], axis=0), reps[0], axis=0)
        return np.concatenate([a, pad], axis=0)
    return a


def normalized_radius(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if r.size == 0:
        return r
    r0 = float(r[0])
    denom = float(r[-1] - r0)
    if abs(denom) < 1e-30:
        return np.zeros_like(r)
    return (r - r0) / denom

# ----------------------------- IMAS helpers -------------------------------

def ensure_imas() -> None:
    if imas is None:
        raise RuntimeError(
            "Could not import imas. Install IMAS-Python / load your IMAS environment. "
            f"Original import error: {IMAS_IMPORT_ERROR!r}"
        )


def make_factory(dd_version: str | None = None):
    ensure_imas()
    if dd_version:
        return imas.IDSFactory(version=dd_version)
    return imas.IDSFactory()


def make_ids(factory: Any, name: str, report: MappingReport):
    if not hasattr(factory, name):
        report.warn(f"Installed IMAS Data Dictionary does not provide IDS '{name}'.")
        return None
    ids = getattr(factory, name)()
    report.bind(ids, name)
    report.created_ids.append(name)
    set_homogeneous_time(ids)
    try:
        ids.ids_properties.comment = (
            "Generated from DREAM HDF5 output by dream_to_imas.py. "
            "Mappings are conservative and reported in the companion mapping report."
        )
    except Exception:
        pass
    return ids


def set_homogeneous_time(ids: Any) -> None:
    try:
        ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    except Exception:
        try:
            ids.ids_properties.homogeneous_time = 1
        except Exception:
            pass


def resize_aos(aos: Any, n: int) -> bool:
    try:
        aos.resize(int(n))
        return True
    except Exception:
        return False


def fill_ids_field(root: Any, path: str, value: Any, report: MappingReport, ids_name: str, source: str) -> bool:
    """Fill one IMAS IDS field and record the mapping in the report.

    Parameters
    ----------
    root:
        IDS object, IDS structure, or AoS item that contains the target field.
    path:
        Slash-separated path relative to ``root``. Numeric components select AoS
        items, for example ``profiles_1d/0/electrons/temperature``.
    value:
        Scalar or array to assign. NumPy scalar types are converted to native
        Python scalars before assignment.
    report:
        Mapping report collecting successful and skipped assignments.
    ids_name:
        Name of the IDS being filled, such as ``"plasma_profiles"``.
    source:
        DREAM source path or derivation text. This is shown in the Markdown
        report, so use precise formulas such as
        ``"volume integral of /eqsys/j_re using current weights"``.

    Returns
    -------
    bool
        ``True`` when the assignment succeeds, ``False`` when the field is absent
        in the installed Data Dictionary or assignment fails.

    This helper is the preferred way to add mapped fields. It keeps the script
    portable across IMAS Data Dictionary versions by recording skipped mappings
    instead of raising on missing nodes.
    """
    parts = [p for p in path.split("/") if p]
    if not parts:
        report.skip(ids_name, path, "empty target path", source, root=root)
        return False
    obj = root
    try:
        for part in parts[:-1]:
            if part.isdigit():
                obj = obj[int(part)]
            else:
                obj = getattr(obj, part)
        leaf = parts[-1]
        if not hasattr(obj, leaf):
            report.skip(ids_name, path, "target node is not present in this DD version", source, root=root)
            return False
        full_node_path = report.full_node_path(ids_name, path, root)
        units = imas_units_from_node(getattr(obj, leaf)) or imas_units_for_node(full_node_path)
        setattr(obj, leaf, sanitize_for_imas(value))
        report.ok(ids_name, path, source, units, root=root)
        return True
    except Exception as exc:
        report.skip(ids_name, path, f"assignment failed: {type(exc).__name__}: {exc}", source, root=root)
        return False


def sanitize_for_imas(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"S", "U", "O"}:
            return value.astype(str)
        return np.asarray(value)
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def resize_child_aos(parent: Any, child_name: str, n: int) -> bool:
    if not hasattr(parent, child_name):
        return False
    return resize_aos(getattr(parent, child_name), n)


# ----------------------------- mappings -----------------------------------

def common_grids(dream: DreamH5, report: MappingReport) -> dict[str, Any]:
    time = flatten_1d(dream.arr("/grid/t", report))

    # Minor radius R-Raxis in [m]
    r  = flatten_1d(dream.arr("/grid/r", report))
    dr = flatten_1d(dream.arr("/grid/dr", report))

    # R at magnetic axis in [m]
    R0 = dream.scalar("/grid/R0", dream.scalar("/grid/eq/R0", None))

    # Flux averaged toroidal field in [T]
    Bphi_avg = flatten_1d(dream.arr("/grid/geometry/GR0", report)) * phi_sign

    # Getting vacuum field from the toroidal field at the outermost radius
    if Bphi_avg is not None and Bphi_avg.size > 0:
        B0 = float(Bphi_avg[-1])
        pass
    
    # Toroidal flux in [Wb]
    phi_tor = flatten_1d(dream.arr("/grid/geometry/toroidalFlux", report))

    rho_tor = np.sqrt(phi_tor / (np.pi * np.abs(B0))) # minor radius like
    rho_tor_norm = normalized_radius(rho_tor) # normalized to 1 at the outermost radius

    # Other geometric quantities
    vpvol  = flatten_1d(dream.arr("/grid/VpVol", report))
    G      = flatten_1d(dream.arr("/grid/geometry/GR0", report)) * phi_sign
    R2inv  = flatten_1d(dream.arr("/grid/geometry/FSA_R02OverR2", report))
    Bmin   = flatten_1d(dream.arr("/grid/geometry/Bmin", report)) * phi_sign
    Bmax   = flatten_1d(dream.arr("/grid/geometry/Bmax", report)) * phi_sign
    FSA_B_over_Bmin = flatten_1d(dream.arr("/grid/geometry/FSA_BOverBmin", report))
    FSA_B_over_Bmin2 = flatten_1d(dream.arr("/grid/geometry/FSA_BOverBmin2", report))
    effective_passing_fraction = flatten_1d(
        dream.arr("/grid/geometry/effectivePassingFraction", report)
    )
    xi0_trapped_boundary = flatten_1d(
        dream.arr("/grid/geometry/xi0TrappedBoundary", report)
    )

    weight_int_area = dr * vpvol * G * R2inv / (2*np.pi*Bmin)
    weight_int_vol  = dr * vpvol * R0
    dVdr =  vpvol * R0

    if time is None or time.size == 0:
        raise ValueError("DREAM file does not contain a valid /grid/t array")
    if r is None or r.size == 0:
        raise ValueError("DREAM file does not contain a valid /grid/r array")
    if dr is None or dr.size == 0:
        raise ValueError("DREAM file does not contain a valid /grid/dr array")

    R_outboard = R0 + r

    return {
        "time": np.asarray(time, dtype=float),
        "r": np.asarray(r, dtype=float),
        "dr": np.asarray(dr, dtype=float),
        "R_outboard": np.asarray(R_outboard, dtype=float),
        "phi_tor": phi_tor,  
        "rho_tor": rho_tor,
        "rho_tor_norm": rho_tor_norm,
        "R0": R0,
        "a": dream.scalar("/grid/a", None),
        "B0": B0,
        "dVdr": dVdr,
        "R2inv": R2inv,
        "Bmin": Bmin,
        "Bmax": Bmax,
        "FSA_B_over_Bmin": FSA_B_over_Bmin,
        "FSA_B_over_Bmin2": FSA_B_over_Bmin2,
        "effective_passing_fraction": effective_passing_fraction,
        "xi0_trapped_boundary": xi0_trapped_boundary,
        "weight_int_area": weight_int_area,
        "weight_int_vol": weight_int_vol
    }


def fill_vacuum_toroidal_field(ids: Any, ids_name: str, grids: dict[str, Any], report: MappingReport) -> None:
    # DREAM output structure shown here does not include B0 explicitly. We fill R0
    # where possible and leave b0 empty unless the user supplies a post-processing hook.
    R0 = grids.get("R0")
    B0 = grids.get("B0")
    if R0 is not None:
        fill_ids_field(ids, "vacuum_toroidal_field/r0", R0, report, ids_name, "/grid/R0")
    else:
        report.warn(
            f"{ids_name}: vacuum_toroidal_field/r0 was not filled because no explicit R0 "
            "dataset was present in the provided DREAM HDF5 structure."
        )
    if B0 is not None:
        time = grids.get("time")
        if time is not None:
            b0 = np.full(np.asarray(time, dtype=float).shape, float(B0), dtype=float)
        else:
            b0 = float(B0)
        fill_ids_field(ids, "vacuum_toroidal_field/b0", b0, report, ids_name, "/grid/B0")
    else:
        report.warn(
            f"{ids_name}: vacuum_toroidal_field/b0 was not filled because no explicit B0 "
            "dataset was present in the provided DREAM HDF5 structure."
        )


def fill_1d_grid(
    p: Any,
    grids: dict[str, Any],
    dream: DreamH5,
    nt: int,
    it: int,
    report: MappingReport,
    ids_name: str,
    psi_p_aligned: Optional[np.ndarray] = None,
) -> None:
    """Fill grid data for a single time step in profiles_1d."""
    rho_tor = grids["rho_tor"]
    rho_tor_norm = grids["rho_tor_norm"]
    fill_ids_field(p, "grid/rho_tor", rho_tor, report, ids_name, "/grid/geometry/toroidalFlux")
    fill_ids_field(p, "grid/rho_tor_norm", rho_tor_norm, report, ids_name, "derived from rho_tor")

    psi_p = psi_p_aligned
    if psi_p is None:
        psi_p = scale_optional(time_aligned(dream.arr("/eqsys/psi_p"), nt), psi_cocos)
    if psi_p is not None:
        psi_arr  = np.asarray(psi_p[it], dtype=float)
        psi_bnd  = float(psi_arr[-1]) if psi_arr.size > 0 else 0.0
        psi_axis = float(psi_arr[0]) if psi_arr.size > 0 else 0.0
        psi_norm = normalized_radius(psi_arr) if psi_arr is not None else None

        fill_ids_field(p, "grid/psi", psi_arr, report, ids_name, "/eqsys/psi_p")
        fill_ids_field(p, "grid/psi_magnetic_axis", psi_axis, report, ids_name, "/eqsys/psi_p[:,0]")
        fill_ids_field(p, "grid/psi_boundary", psi_bnd, report, ids_name, "/eqsys/psi_p[:,-1]")
        
        fill_ids_field(p, "grid/rho_pol_norm", psi_norm, report, ids_name, "derived")


def map_plasma_profiles(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "plasma_profiles"
    pp = make_ids(factory, ids_name, report)
    if pp is None:
        return None

    time = grids["time"]
    nt = len(time)

    fill_ids_field(pp, "time", time, report, ids_name, "/grid/t")
    fill_vacuum_toroidal_field(pp, ids_name, grids, report)

    if not resize_child_aos(pp, "profiles_1d", nt):
        report.skip(ids_name, "profiles_1d", "could not resize AoS", "")
        return pp

    # DREAM time-dependent profiles.
    profiles = {
        "electrons/temperature": ("/eqsys/T_cold", dream.arr("/eqsys/T_cold", report)),
        "electrons/density": ("/eqsys/n_tot", dream.arr("/eqsys/n_tot")),
        "electrons/density_thermal": ("/eqsys/n_cold", dream.arr("/eqsys/n_cold", report)),
        "conductivity_parallel": ("/other/fluid/conductivity", dream.arr("/other/fluid/conductivity")),
        "e_field/parallel": ("/eqsys/E_field", dream.arr("/eqsys/E_field", report)),
        "zeff": ("/other/fluid/Zeff", dream.arr("/other/fluid/Zeff")),
    }

    # n_tot is total electron density; if absent, use n_cold as the main density.
    if profiles["electrons/density"][1] is None:
        profiles["electrons/density"] = ("/eqsys/n_cold", dream.arr("/eqsys/n_cold", report))

    j_tot = dream.arr("/eqsys/j_tot", report) * phi_sign
    j_ohm = dream.arr("/eqsys/j_ohm", report) * phi_sign
    aligned_profiles = {target: (source, time_aligned(data, nt)) for target, (source, data) in profiles.items()}
    j_tot_aligned = time_aligned(j_tot, nt)
    j_ohm_aligned = time_aligned(j_ohm, nt)
    psi_p_aligned = scale_optional(time_aligned(dream.arr("/eqsys/psi_p"), nt), psi_cocos)
    ion_context = prepare_ion_profile_context(dream, nt, report)

    for it in range(nt):
        p = pp.profiles_1d[it]
        report.bind(p, "plasma_profiles.profiles_1d")
        fill_ids_field(p, "time", time[it], report, ids_name, "/grid/t")

        fill_1d_grid(p, grids, dream, nt, it, report, ids_name, psi_p_aligned)

        for target, (source, aligned) in aligned_profiles.items():
            if aligned is not None:
                fill_ids_field(p, target, aligned[it], report, ids_name, source)
            elif source:
                report.missing(source)

        if j_tot_aligned is not None:
            fill_ids_field(p, "j_total", j_tot_aligned[it], report, ids_name, "/eqsys/j_tot")
        if j_ohm_aligned is not None:
            fill_ids_field(p, "j_ohmic", j_ohm_aligned[it], report, ids_name, "/eqsys/j_ohm")

        # Ion densities: DREAM /eqsys/n_i has shape (time, charge_state_flat, r).
        fill_ion_profiles_1d(p, dream, it, nt, report, ids_name, ion_context)

    return pp


def dream_arr_first(dream: Any, *paths: str) -> Any:
    """
    Return the first available DREAM HDF5 dataset among `paths`.

    This lets us prefer /settings/eqsys/n_i/*, while falling back to /ionmeta/*.
    """
    for path in paths:
        try:
            value = dream.arr(path)
        except Exception:
            value = None
        if value is not None:
            return value
    return None


def parse_dream_string_list(value: Any) -> list[str]:
    """
    Parse DREAM string-list datasets.

    DREAM often stores lists as one character array, for example:
        D;T;D_inj_stage_1;Ne_inj_stage_1;

    This returns:
        ["D", "T", "D_inj_stage_1", "Ne_inj_stage_1"]
    """
    if value is None:
        return []

    if isinstance(value, str):
        text = value

    elif isinstance(value, (bytes, bytearray)):
        text = bytes(value).decode("utf-8", errors="ignore")

    else:
        arr = np.asarray(value)

        if arr.dtype.kind in {"S", "U", "O"}:
            chars = []
            for x in arr.ravel():
                if isinstance(x, bytes):
                    chars.append(x.decode("utf-8", errors="ignore"))
                else:
                    chars.append(str(x))
            text = "".join(chars)

        elif arr.dtype.kind in {"u", "i"}:
            text = bytes(arr.ravel().astype(np.uint8)).decode("utf-8", errors="ignore")

        else:
            text = str(value)

    text = text.replace("\x00", "").strip()

    # DREAM convention: semicolon-separated list.
    if ";" in text:
        return [p.strip() for p in text.split(";") if p.strip()]

    # Fallback for already-separated strings.
    for sep in [",", "\n"]:
        if sep in text:
            return [p.strip() for p in text.split(sep) if p.strip()]

    return [text] if text else []


def parse_dream_ion_names(value: Any, n: int | None = None) -> list[str]:
    """
    Backwards-compatible wrapper around parse_dream_string_list().
    """
    names = parse_dream_string_list(value)
    if n is not None and len(names) > n:
        # Do not blindly truncate here. The reconciliation routine below will decide
        # which names are compatible with Z/A.
        return names
    return names


def base_species_label(label: str) -> str:
    """
    Extract the physical species prefix from DREAM labels.

    Examples:
        D                -> D
        T                -> T
        D_inj_stage_1    -> D
        Ne_inj_stage_2   -> Ne
    """
    head = label.split("_", 1)[0]
    if len(head) >= 2 and head[:2].istitle():
        return head[:2]
    return head[:1]


def infer_isotope_mass_number(
    label: str,
    z: int,
    isotope: int | None,
    hydrogen_names: set[str],
    tritium_names: set[str],
) -> int | None:
    """
    Infer isotope mass number A.

    Priority:
      1. /settings/eqsys/n_i/isotopes, if present and positive
      2. explicit tritiumnames / hydrogennames
      3. label prefix H/D/T
      4. DREAM convention: Z=1 defaults to deuterium
    """
    if isotope is not None and isotope > 0:
        return int(isotope)

    if z == 1:
        if label in hydrogen_names:
            return 1
        if label in tritium_names:
            return 3

        base = base_species_label(label)
        if base == "H":
            return 1
        if base == "D":
            return 2
        if base == "T":
            return 3

        # DREAM default for Z=1 if not explicitly marked otherwise.
        return 2
    
    if z == 10:
        return 20   # Neon
    if z == 18:
        return 40   # Argon

    return None


def element_label_from_z(z: int) -> str:
    element = {
        1: "H",
        2: "He",
        3: "Li",
        4: "Be",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "Ne",
        18: "Ar",
        74: "W",
    }
    return element.get(int(z), f"Z{int(z)}")


def expected_z_from_label(label: str) -> int | None:
    """
    Return nuclear charge inferred from the label prefix, when obvious.

    This is used only to avoid assigning e.g. a D label to a Ne block.
    """
    element_z = {
        "H": 1,
        "D": 1,
        "T": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Ar": 18,
        "W": 74,
    }

    base = base_species_label(label)
    return element_z.get(base)


def reconcile_dream_ion_names(
    raw_names: Any,
    z_list: list[int],
    isotope_list: list[int | None],
    report: Any | None = None,
) -> list[str]:
    """
    Return exactly len(z_list) labels.

    /settings/eqsys/n_i/Z and /eqsys/n_i define the real species blocks.
    /settings/eqsys/n_i/names or /ionmeta/names provide labels, but in some
    injection cases there can be extra labels. We therefore only accept labels
    compatible with the species Z.
    """
    parsed = parse_dream_ion_names(raw_names)
    out: list[str] = []
    used: set[int] = set()

    for i, z in enumerate(z_list):
        chosen = None

        # First try same-position label.
        if i < len(parsed):
            z_label = expected_z_from_label(parsed[i])
            if z_label is None or z_label == z:
                chosen = parsed[i]
                used.add(i)

        # If same-position failed, scan unused labels for a compatible one.
        if chosen is None:
            for j, label in enumerate(parsed):
                if j in used:
                    continue
                z_label = expected_z_from_label(label)
                if z_label is None or z_label == z:
                    chosen = label
                    used.add(j)
                    break

        if chosen is None:
            A = isotope_list[i] if i < len(isotope_list) else None
            chosen = f"species_{i}_Z{z}" if A is None else f"species_{i}_Z{z}_A{A}"
            if report is not None:
                report.warn(
                    f"Could not find a DREAM ion name compatible with species {i}, "
                    f"Z={z}. Parsed names were {parsed}. Using '{chosen}'."
                )

        out.append(chosen)

    if report is not None and len(parsed) != len(z_list):
        report.warn(
            f"Parsed {len(parsed)} DREAM ion names but found {len(z_list)} ion species "
            f"from Z/isotopes. Names were reconciled by compatibility with Z."
        )

    if report is not None:
        unused = [label for j, label in enumerate(parsed) if j not in used]
        if unused:
            report.warn(f"Unused DREAM ion labels not mapped to n_i species blocks: {unused}")

    return out


def set_element_atomic_properties(
    parent: Any,
    label: str,
    z: int,
    A: int | None,
    report: Any,
    ids_name: str,
    source: str,
) -> None:
    """
    Store atomic properties in parent.element[0], where available.

    Intended for both:
        profiles_1d.ion[i].element[0]
        profiles_1d.neutral[i].element[0]
    """
    if not hasattr(parent, "element") or not resize_aos(parent.element, 1):
        report.skip(ids_name, "element", "element AoS not available/resizable", source, root=parent)
        return

    element = parent.element[0]
    parent_path = report.object_paths.get(id(parent))
    if parent_path:
        report.bind(element, f"{parent_path}.element")

    # IMAS commonly uses z_n for nuclear charge and a for mass number.
    fill_ids_field(element, "z_n", float(z), report, ids_name, source)

    if A is not None:
        fill_ids_field(element, "a", float(A), report, ids_name, source)


@dataclass
class IonProfileContext:
    n_i: np.ndarray
    z_list: list[int]
    labels: list[str]
    a_list: list[int | None]


def prepare_ion_profile_context(dream: Any, nt: int, report: Any) -> Optional[IonProfileContext]:
    n_i = time_aligned(dream_arr_first(dream, "/eqsys/n_i"), nt)
    z_raw = dream_arr_first(dream, "/settings/eqsys/n_i/Z", "/ionmeta/Z")
    if n_i is None or z_raw is None:
        return None

    z_arr = flatten_1d(z_raw)
    if z_arr is None:
        return None
    z_list = [int(z) for z in z_arr]

    isotopes_raw = dream_arr_first(dream, "/settings/eqsys/n_i/isotopes")
    isotope_arr = flatten_1d(isotopes_raw) if isotopes_raw is not None else None
    isotope_list: list[int | None] = []
    for i in range(len(z_list)):
        if isotope_arr is not None and i < len(isotope_arr):
            mass = int(isotope_arr[i])
            isotope_list.append(mass if mass > 0 else None)
        else:
            isotope_list.append(None)

    names_raw = dream_arr_first(dream, "/settings/eqsys/n_i/names", "/ionmeta/names")
    labels = reconcile_dream_ion_names(
        names_raw,
        z_list=z_list,
        isotope_list=isotope_list,
        report=report,
    )

    hydrogen_names = set(parse_dream_string_list(dream_arr_first(dream, "/settings/eqsys/n_i/hydrogennames")))
    tritium_names = set(parse_dream_string_list(dream_arr_first(dream, "/settings/eqsys/n_i/tritiumnames")))

    a_list = [
        infer_isotope_mass_number(
            label=labels[i],
            z=z_list[i],
            isotope=isotope_list[i],
            hydrogen_names=hydrogen_names,
            tritium_names=tritium_names,
        )
        for i in range(len(z_list))
    ]

    expected = sum(z + 1 for z in z_list)
    if n_i.shape[1] < expected:
        report.warn(
            f"/eqsys/n_i has {n_i.shape[1]} charge-state channels, but "
            f"Z implies {expected}; mapping the available prefix only."
        )
    elif n_i.shape[1] > expected:
        report.warn(
            f"/eqsys/n_i has {n_i.shape[1]} charge-state channels, but "
            f"Z implies {expected}; extra trailing channels will not be mapped."
        )

    return IonProfileContext(n_i=n_i, z_list=z_list, labels=labels, a_list=a_list)


def fill_ion_profiles_1d(
    p: Any,
    dream: Any,
    it: int,
    nt: int,
    report: Any,
    ids_name: str,
    context: Optional[IonProfileContext] = None,
) -> None:
    """
    Map DREAM ion densities to IMAS plasma_profiles.profiles_1d.

    DREAM:
        /eqsys/n_i has shape (time, charge_state_flat, radius)

    DREAM ordering:
        species 0, Z0=0
        species 0, Z0=1
        ...
        species 0, Z0=Z

        species 1, Z0=0
        species 1, Z0=1
        ...

    Mapping:
        Z0=0      -> p.neutral[i].density
        Z0=1..Z   -> p.ion[i].state[Z0-1].density
        ion.density = sum over charged states only

    Atomic properties:
        p.ion[i].element[0].z_n     = Z
        p.ion[i].element[0].a       = A
        p.neutral[i].element[0].z_n = Z
        p.neutral[i].element[0].a   = A
    """
    context = context or prepare_ion_profile_context(dream, nt, report)
    if context is None:
        return

    n_i = context.n_i
    z_list = context.z_list
    labels = context.labels
    A_list = context.a_list

    has_ion = hasattr(p, "ion") and resize_aos(p.ion, len(z_list))
    if not has_ion:
        report.skip(ids_name, "profiles_1d/ion", "ion AoS not available/resizable", "/eqsys/n_i")

    has_neutral = hasattr(p, "neutral") and resize_aos(p.neutral, len(z_list))
    if not has_neutral:
        report.skip(ids_name, "profiles_1d/neutral", "neutral AoS not available/resizable", "/eqsys/n_i")

    if not has_ion and not has_neutral:
        return

    idx = 0

    for iion, z in enumerate(z_list):
        label = labels[iion]
        A = A_list[iion]
        nstates = z + 1

        block = n_i[it, idx:min(idx + nstates, n_i.shape[1]), :]
        idx += nstates

        if block.size == 0:
            continue

        if has_ion and has_neutral:
            report.bind(p.ion[iion], "plasma_profiles.profiles_1d.ion")
            report.bind(p.neutral[iion], "plasma_profiles.profiles_1d.neutral")
            fill_ids_field(p.ion[iion], "neutral_index", iion, report, ids_name, "DREAM species index")
            fill_ids_field(p.neutral[iion], "ion_index", iion, report, ids_name, "DREAM species index")

        # ------------------------------------------------------------
        # DREAM Z0=0: neutral density
        # ------------------------------------------------------------
        neutral_density = block[0, :]

        if has_neutral:
            neutral = p.neutral[iion]
            report.bind(neutral, "plasma_profiles.profiles_1d.neutral")

            set_element_atomic_properties(
                neutral,
                label=label,
                z=z,
                A=A,
                report=report,
                ids_name=ids_name,
                source="/settings/eqsys/n_i/Z,/settings/eqsys/n_i/isotopes",
            )

            # Depending on DD version, neutral density may be direct or under state.
            fill_ids_field(neutral, "density", neutral_density, report, ids_name, "/eqsys/n_i Z0=0")

            if hasattr(neutral, "state") and resize_aos(neutral.state, 1):
                nstate = neutral.state[0]
                report.bind(nstate, "plasma_profiles.profiles_1d.neutral.state")
                fill_ids_field(nstate, "density", neutral_density, report, ids_name, "/eqsys/n_i Z0=0")

        # ------------------------------------------------------------
        # DREAM Z0=1..Z: charged ion densities
        # ------------------------------------------------------------
        charged_block = block[1:, :]

        if has_ion:
            ion = p.ion[iion]
            report.bind(ion, "plasma_profiles.profiles_1d.ion")

            fill_ids_field(ion, "name", label, report, ids_name, "/settings/eqsys/n_i/names")
            set_element_atomic_properties(
                ion,
                label=label,
                z=z,
                A=A,
                report=report,
                ids_name=ids_name,
                source="/settings/eqsys/n_i/Z,/settings/eqsys/n_i/isotopes",
            )

            # ion.density should exclude neutrals.
            if charged_block.size:
                fill_ids_field(
                    ion,
                    "density",
                    np.sum(charged_block, axis=0),
                    report,
                    ids_name,
                    "/eqsys/n_i summed over charged states Z0=1..Z",
                )

            # One IMAS ion.state entry per charged state.
            # DREAM block[1] -> Z0=1, block[2] -> Z0=2, ...
            if hasattr(ion, "state") and resize_aos(ion.state, charged_block.shape[0]):
                for z0 in range(1, block.shape[0]):
                    state = ion.state[z0 - 1]
                    report.bind(state, "plasma_profiles.profiles_1d.ion.state")

                    # Do NOT set z_min/z_max here, per your request.
                    # The charge state is encoded by ordering and label.
                    fill_ids_field(
                        state,
                        "name",
                        f"{label} Z0={z0}",
                        report,
                        ids_name,
                        "/eqsys/n_i charge-state index",
                    )
                    fill_ids_field(
                        state,
                        "density",
                        block[z0, :],
                        report,
                        ids_name,
                        "/eqsys/n_i",
                    )

def map_runaway_electrons(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "runaway_electrons"
    re_ids = make_ids(factory, ids_name, report)
    if re_ids is None:
        return None

    time = grids["time"]
    nt   = len(time)

    fill_ids_field(re_ids, "time", time, report, ids_name, "/grid/t")
    fill_vacuum_toroidal_field(re_ids, ids_name, grids, report)

    if not resize_child_aos(re_ids, "profiles_1d", nt):
        report.skip(ids_name, "profiles_1d", "could not resize AoS", "")
        return re_ids

    quantities = {
        "density": ("/eqsys/n_re", dream.arr("/eqsys/n_re", report)),
        "current_density": ("/eqsys/j_re", dream.arr("/eqsys/j_re", report) * phi_sign),
        #"energy_density_kinetic": ("/other/fluid/W_re", dream.arr("/other/fluid/W_re")),
        "ddensity_dt_total": ("/other/fluid/runawayRate", dream.arr("/other/fluid/runawayRate")),
        "ddensity_dt_dreicer": ("/other/fluid/gammaDreicer", dream.arr("/other/fluid/gammaDreicer")),
        "ddensity_dt_hot_tail": ("/other/fluid/gammaHottail", dream.arr("/other/fluid/gammaHottail")), 
        "ddensity_dt_tritium": ("/other/fluid/gammaTritium", dream.arr("/other/fluid/gammaTritium")),  
        "ddensity_dt_compton": ("/other/fluid/gammaCompton", dream.arr("/other/fluid/gammaCompton")),  
        #"ddensity_dt_avalanche": ("/other/fluid/GammaAva", dream.arr("/other/fluid/GammaAva")),  #not available in dictionary
        "momentum_critical_avalanche": ("/other/fluid/pCrit", scale_optional(dream.arr("/other/fluid/pCrit"), p_norm)),
        "momentum_critical_hot_tail": ("/other/fluid/pCritHottail", scale_optional(dream.arr("/other/fluid/pCritHottail"), p_norm)),
        "e_field_dreicer": ("/other/fluid/EDreic", dream.arr("/other/fluid/EDreic")),
        "e_field_critical": ("/other/fluid/Ectot", dream.arr("/other/fluid/Ectot")),
    }
    aligned_quantities = {
        target: (source, time_aligned(data, nt))
        for target, (source, data) in quantities.items()
    }
    psi_p_aligned = scale_optional(time_aligned(dream.arr("/eqsys/psi_p"), nt), psi_cocos)

    # Fill up profiles_1d for each time step
    for it in range(nt):
        p = re_ids.profiles_1d[it]
        report.bind(p, "runaway_electrons.profiles_1d")
        fill_ids_field(p, "time", time[it], report, ids_name, "/grid/t")

        fill_1d_grid(p, grids, dream, nt, it, report, ids_name, psi_p_aligned)

        for target, (source, aligned) in aligned_quantities.items():
            if aligned is not None:
                fill_ids_field(p, target, aligned[it], report, ids_name, source)
            elif source:
                report.missing(source)


    j_re = dream.arr("/eqsys/j_re") * phi_sign
    I_re = current_from_j_trace(j_re, grids, nt)

    fill_ids_field(re_ids, "global_quantities/current_phi", I_re, report, ids_name, "derived from j_re")

    return re_ids


def reshape_spi_vector(data: Optional[np.ndarray], nt: int) -> Optional[np.ndarray]:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    if arr.ndim == 2 and arr.shape[1] % 3 == 0:
        return arr.reshape((arr.shape[0], arr.shape[1] // 3, 3))
    return None


def reshape_spi_scalar(data: Optional[np.ndarray], nt: int) -> Optional[np.ndarray]:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    if arr.ndim == 2:
        return arr
    return None


def dream_spi_position_to_imas_rphiz(xyz: np.ndarray, R0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert DREAM local Cartesian SPI positions to IMAS cylindrical coordinates.

    DREAM stores fragment positions as local Cartesian coordinates relative to
    the magnetic axis. The first component is the major-radius offset, the
    second component is vertical, and the third component is the Cartesian
    toroidal displacement. IMAS SPI expects cylindrical R/phi/Z.
    """
    xyz = np.asarray(xyz, dtype=float)
    x = xyz[..., 0]
    z = xyz[..., 1]
    y_tor = xyz[..., 2] if xyz.shape[-1] > 2 else np.zeros_like(x)
    r_cart = R0 + x
    r = np.hypot(r_cart, y_tor)
    phi = np.arctan2(y_tor, r_cart)
    return r, phi, z


def dream_spi_velocity_to_imas_rphiz(
    xyz: np.ndarray,
    vxyz: np.ndarray,
    R0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert DREAM local Cartesian SPI velocities to cylindrical components."""
    xyz = np.asarray(xyz, dtype=float)
    vxyz = np.asarray(vxyz, dtype=float)
    x = xyz[..., 0]
    y_tor = xyz[..., 2] if xyz.shape[-1] > 2 else np.zeros_like(x)
    vx = vxyz[..., 0]
    vz = vxyz[..., 1]
    vy_tor = vxyz[..., 2] if vxyz.shape[-1] > 2 else np.zeros_like(vx)
    r_cart = R0 + x
    r = np.hypot(r_cart, y_tor)
    r_safe = np.where(r > 0.0, r, 1.0)
    velocity_r = (r_cart * vx + y_tor * vy_tor) / r_safe
    velocity_phi = (r_cart * vy_tor - y_tor * vx) / r_safe
    return velocity_r, velocity_phi, vz


@dataclass
class SpiComponent:
    label: str
    z: int
    isotope: int | None
    atoms_n: float | None
    source: str


@dataclass
class SpiInjectorGroup:
    shard_indices: np.ndarray
    fractions: np.ndarray | None
    source: str


def spi_display_label(label: str, z: int, isotope: int | None) -> str:
    base = base_species_label(label) or element_label_from_z(z)
    if z == 1 and isotope == 1:
        return "H"
    if z == 1 and isotope == 2:
        return "D"
    if z == 1 and isotope == 3:
        return "T"
    return base if expected_z_from_label(base) == z else element_label_from_z(z)


def spi_components_from_settings(dream: DreamH5, report: MappingReport) -> list[SpiComponent]:
    atoms_raw = dream.arr("/settings/eqsys/spi/init/Ninj")
    atoms_arr = flatten_1d(atoms_raw)
    n_components = int(atoms_arr.size) if atoms_arr is not None else 0

    names_raw = dream.arr("/settings/eqsys/n_i/names")
    names = parse_dream_string_list(names_raw)
    injected_names = [(i, name) for i, name in enumerate(names) if "_inj" in name]

    n_i_z = flatten_1d(dream.arr("/settings/eqsys/n_i/Z"))
    n_i_isotopes = flatten_1d(dream.arr("/settings/eqsys/n_i/isotopes"))

    spi_z = flatten_1d(dream.arr("/settings/eqsys/spi/ZsDrift"))
    spi_isotopes = flatten_1d(dream.arr("/settings/eqsys/spi/isotopesDrift"))

    if n_components == 0 and spi_z is not None:
        n_components = int(spi_z.size)
        atoms_arr = np.full((n_components,), np.nan)

    components: list[SpiComponent] = []
    use_injected_names = n_components > 0 and len(injected_names) == n_components

    for i in range(n_components):
        atoms = None
        if atoms_arr is not None and i < atoms_arr.size and np.isfinite(atoms_arr[i]):
            atoms = float(atoms_arr[i])

        label = ""
        z: int | None = None
        isotope: int | None = None
        source = "/settings/eqsys/spi/init/Ninj"

        if use_injected_names:
            name_index, label = injected_names[i]
            source = "/settings/eqsys/n_i/names, /settings/eqsys/spi/init/Ninj"
            if n_i_z is not None and name_index < n_i_z.size:
                z = int(n_i_z[name_index])
            if n_i_isotopes is not None and name_index < n_i_isotopes.size:
                raw_isotope = int(n_i_isotopes[name_index])
                isotope = raw_isotope if raw_isotope > 0 else None

        if z is None and label:
            z = expected_z_from_label(label)

        if z is None and spi_z is not None and spi_z.size > 0:
            z = int(spi_z[min(i, spi_z.size - 1)])
            source = "/settings/eqsys/spi/ZsDrift, /settings/eqsys/spi/init/Ninj"

        if z is None:
            report.warn(f"Could not infer SPI component {i} species; skipping it.")
            continue

        if isotope is None and spi_z is not None and spi_isotopes is not None:
            matches = np.flatnonzero(np.asarray(spi_z, dtype=int) == z)
            if matches.size > 0 and matches[0] < spi_isotopes.size:
                raw_isotope = int(spi_isotopes[matches[0]])
                isotope = raw_isotope if raw_isotope > 0 else None

        if not label:
            label = element_label_from_z(z)

        isotope = infer_isotope_mass_number(label, z, isotope, set(), set())
        components.append(SpiComponent(label=label, z=z, isotope=isotope, atoms_n=atoms, source=source))

    return components


def spi_molar_fraction_matrix(
    dream: DreamH5,
    n_components: int,
    n_shards: int,
    report: MappingReport,
) -> np.ndarray | None:
    raw = flatten_1d(dream.arr("/settings/eqsys/n_i/SPIMolarFraction"))
    if raw is None or n_components <= 0 or n_shards <= 0:
        return None

    expected = n_components * n_shards
    if raw.size == expected + 1 and raw[0] < 0:
        values = raw[1:]
    elif raw.size == expected:
        values = raw
    else:
        report.warn(
            "/settings/eqsys/n_i/SPIMolarFraction has incompatible size "
            f"{raw.size}; expected {expected} or {expected + 1}. Falling back to one SPI injector."
        )
        return None

    return np.asarray(values, dtype=float).reshape(n_components, n_shards)


def spi_injector_groups(dream: DreamH5, components: list[SpiComponent], n_shards: int, report: MappingReport) -> list[SpiInjectorGroup]:
    fractions = spi_molar_fraction_matrix(dream, len(components), n_shards, report)
    if fractions is None:
        return [
            SpiInjectorGroup(
                shard_indices=np.arange(n_shards, dtype=int),
                fractions=None,
                source="/settings/eqsys/spi/ZsDrift, /settings/eqsys/spi/init/Ninj",
            )
        ]

    shard_fractions = fractions.T
    rounded = np.round(shard_fractions, 9)
    groups: list[SpiInjectorGroup] = []
    start = 0
    for ish in range(1, n_shards):
        if not np.array_equal(rounded[ish], rounded[ish - 1]):
            group_fraction = np.mean(shard_fractions[start:ish], axis=0)
            groups.append(
                SpiInjectorGroup(
                    shard_indices=np.arange(start, ish, dtype=int),
                    fractions=group_fraction,
                    source="/settings/eqsys/n_i/SPIMolarFraction",
                )
            )
            start = ish

    group_fraction = np.mean(shard_fractions[start:n_shards], axis=0)
    groups.append(
        SpiInjectorGroup(
            shard_indices=np.arange(start, n_shards, dtype=int),
            fractions=group_fraction,
            source="/settings/eqsys/n_i/SPIMolarFraction",
        )
    )

    valid_groups = [group for group in groups if group.shard_indices.size > 0 and np.any(np.abs(group.fractions) > 1e-12)]
    if not valid_groups:
        report.warn("SPI molar fraction data did not contain non-zero shard compositions; falling back to one injector.")
        return [
            SpiInjectorGroup(
                shard_indices=np.arange(n_shards, dtype=int),
                fractions=None,
                source="/settings/eqsys/spi/ZsDrift, /settings/eqsys/spi/init/Ninj",
            )
        ]
    return valid_groups


def spi_group_component_weights(
    groups: list[SpiInjectorGroup],
    n_components: int,
    initial_volume: np.ndarray | None,
) -> np.ndarray:
    weights = np.zeros((len(groups), n_components), dtype=float)
    for igroup, group in enumerate(groups):
        if group.fractions is None:
            weights[igroup, :] = 1.0
            continue
        shard_weight = np.ones(group.shard_indices.size, dtype=float)
        if initial_volume is not None and initial_volume.size > int(np.max(group.shard_indices)):
            shard_weight = np.asarray(initial_volume[group.shard_indices], dtype=float)
        for icomp in range(n_components):
            weights[igroup, icomp] = float(np.sum(shard_weight * max(float(group.fractions[icomp]), 0.0)))
    return weights


def map_spi_species(
    parent: Any,
    components: list[SpiComponent],
    group: SpiInjectorGroup,
    component_weights: np.ndarray,
    group_index: int,
    initial_volume: np.ndarray | None,
    report: MappingReport,
    ids_name: str,
    source_prefix: str,
) -> None:
    if not components:
        return

    active: list[tuple[int, SpiComponent, float, float | None]] = []
    if group.fractions is None:
        for icomp, component in enumerate(components):
            atoms = component.atoms_n
            fraction = None
            active.append((icomp, component, 1.0, atoms))
    else:
        component_totals = np.sum(component_weights, axis=0)
        for icomp, component in enumerate(components):
            fraction = float(group.fractions[icomp])
            if abs(fraction) <= 1e-12:
                continue
            atoms = None
            if component.atoms_n is not None and icomp < component_totals.size and component_totals[icomp] > 0.0:
                atoms = float(component.atoms_n * component_weights[group_index, icomp] / component_totals[icomp])
            active.append((icomp, component, fraction, atoms))

    atoms_values = [atoms for _, _, _, atoms in active if atoms is not None]
    total_atoms = float(np.sum(atoms_values)) if atoms_values and np.sum(atoms_values) > 0.0 else None
    if total_atoms is not None:
        fill_ids_field(parent, "atoms_n", total_atoms, report, ids_name, "/settings/eqsys/spi/init/Ninj")

    total_volume = None
    if initial_volume is not None and group.shard_indices.size > 0:
        max_index = int(np.max(group.shard_indices))
        if initial_volume.size > max_index:
            volume_sum = float(np.sum(initial_volume[group.shard_indices]))
            if volume_sum > 0.0:
                total_volume = volume_sum

    if not resize_child_aos(parent, "species", len(active)):
        report.skip(ids_name, f"{source_prefix}/species", "could not resize AoS", group.source)
        return

    fraction_sum = None
    if group.fractions is not None:
        fraction_sum = float(np.sum([max(fraction, 0.0) for _, _, fraction, _ in active]))

    for ispecies, (_, component, fraction, atoms) in enumerate(active):
        label = spi_display_label(component.label, component.z, component.isotope)
        species = parent.species[ispecies]
        parent_path = report.object_paths.get(id(parent), "spi.injector.pellet.core")
        report.bind(species, f"{parent_path}.species")
        fill_ids_field(species, "name", label, report, ids_name, component.source)
        fill_ids_field(species, "z_n", float(component.z), report, ids_name, component.source)
        if component.isotope is not None:
            fill_ids_field(species, "a", float(component.isotope), report, ids_name, component.source)

        mixture_fraction = None
        if group.fractions is not None and fraction_sum is not None and fraction_sum > 0.0:
            mixture_fraction = float(max(fraction, 0.0) / fraction_sum)
        elif atoms is not None and total_atoms is not None:
            mixture_fraction = float(atoms / total_atoms)

        if mixture_fraction is not None and total_atoms is not None and total_volume is not None:
            density = float(total_atoms / total_volume * mixture_fraction)
            fill_ids_field(species, "density", density, report, ids_name, f"{group.source}, /eqsys/Y_p or /settings/eqsys/spi/init/rp")
        elif mixture_fraction is not None:
            report.skip(
                ids_name,
                f"{source_prefix}/species/{ispecies}/density",
                "cannot compute atomic density without pellet volume and atoms_n",
                f"{group.source}, /settings/eqsys/spi/init/Ninj",
            )


def map_spi(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "spi"
    spi = make_ids(factory, ids_name, report)
    if spi is None:
        return None

    time = grids["time"]
    nt = len(time)
    R0 = float(grids.get("R0") or 0.0)
    fill_ids_field(spi, "time", time, report, ids_name, "/grid/t")

    x_p = reshape_spi_vector(dream.arr("/eqsys/x_p", report), nt)
    v_p = reshape_spi_vector(dream.arr("/eqsys/v_p", report), nt)
    y_p = reshape_spi_scalar(dream.arr("/eqsys/Y_p", report), nt)

    if x_p is None:
        x_init = dream.arr("/settings/eqsys/spi/init/xp", report)
        if x_init is not None:
            x_init = np.asarray(x_init, dtype=float).reshape(-1)
            if x_init.size % 3 == 0:
                x_p = np.repeat(x_init.reshape(1, x_init.size // 3, 3), nt, axis=0)
    if v_p is None:
        v_init = dream.arr("/settings/eqsys/spi/init/vp", report)
        if v_init is not None:
            v_init = np.asarray(v_init, dtype=float).reshape(-1)
            if v_init.size % 3 == 0:
                v_p = np.repeat(v_init.reshape(1, v_init.size // 3, 3), nt, axis=0)
    if y_p is None:
        rp_init = dream.arr("/settings/eqsys/spi/init/rp", report)
        if rp_init is not None:
            rp_init = np.asarray(rp_init, dtype=float).reshape(-1)
            y_p = np.repeat((rp_init**(5.0 / 3.0)).reshape(1, -1), nt, axis=0)

    if x_p is None and y_p is None:
        report.skip(ids_name, "injector", "no DREAM SPI shard position or radius data found", "/eqsys/x_p,/eqsys/Y_p")
        return spi

    nshards = x_p.shape[1] if x_p is not None else y_p.shape[1]

    volume = None
    if y_p is not None:
        radius = np.maximum(y_p, 0.0) ** (3.0 / 5.0)
        volume = (4.0 / 3.0) * np.pi * radius**3

    components = spi_components_from_settings(dream, report)
    groups = spi_injector_groups(dream, components, nshards, report)
    initial_volume = volume[0] if volume is not None and volume.shape[0] > 0 else None
    component_weights = spi_group_component_weights(groups, len(components), initial_volume)

    if not resize_child_aos(spi, "injector", len(groups)):
        report.skip(ids_name, "injector", "could not resize AoS", "")
        return spi

    for igroup, group in enumerate(groups):
        injector = spi.injector[igroup]
        report.bind(injector, "spi.injector")
        name = "DREAM_SPI" if len(groups) == 1 else f"DREAM_SPI_{igroup + 1}"
        fill_ids_field(injector, "name", name, report, ids_name, "DREAM SPI")
        fill_ids_field(
            injector,
            "description",
            "Shattered pellet injection reconstructed from DREAM shard state",
            report,
            ids_name,
            "DREAM SPI",
        )

        if hasattr(injector, "pellet") and hasattr(injector.pellet, "core"):
            report.bind(injector.pellet.core, "spi.injector.pellet.core")
            map_spi_species(
                injector.pellet.core,
                components,
                group,
                component_weights,
                igroup,
                initial_volume,
                report,
                ids_name,
                f"injector/{igroup}/pellet/core",
            )

        group_shards = group.shard_indices

        if x_p is not None and x_p.shape[1] > int(np.max(group_shards)):
            r0, phi0, z0 = dream_spi_position_to_imas_rphiz(x_p[0, group_shards, :], R0)
            fill_ids_field(injector, "shattering_position/r", float(np.mean(r0)), report, ids_name, "/eqsys/x_p[0]")
            fill_ids_field(injector, "shattering_position/phi", float(np.mean(phi0)), report, ids_name, "/eqsys/x_p[0]")
            fill_ids_field(injector, "shattering_position/z", float(np.mean(z0)), report, ids_name, "/eqsys/x_p[0]")
            fill_ids_field(injector, "shatter_cone/origin/r", float(np.mean(r0)), report, ids_name, "/eqsys/x_p[0]")
            fill_ids_field(injector, "shatter_cone/origin/phi", float(np.mean(phi0)), report, ids_name, "/eqsys/x_p[0]")
            fill_ids_field(injector, "shatter_cone/origin/z", float(np.mean(z0)), report, ids_name, "/eqsys/x_p[0]")

        # this should be set only for the very first DREAM file of the simulation
        if (
            x_p is not None
            and v_p is not None
            and x_p.shape[1] > int(np.max(group_shards))
            and v_p.shape[1] > int(np.max(group_shards))
        ):
            vr0, vphi0, vz0 = dream_spi_velocity_to_imas_rphiz(x_p[0, group_shards, :], v_p[0, group_shards, :], R0)
            fill_ids_field(
                injector,
                "velocity_mass_centre_fragments_r",
                float(np.mean(vr0)),
                report,
                ids_name,
                "cylindrical transform of /eqsys/x_p[0] and /eqsys/v_p[0]",
            )
            fill_ids_field(
                injector,
                "velocity_mass_centre_fragments_phi",
                float(np.mean(vphi0)),
                report,
                ids_name,
                "cylindrical transform of /eqsys/x_p[0] and /eqsys/v_p[0]",
            )
            fill_ids_field(
                injector,
                "velocity_mass_centre_fragments_z",
                float(np.mean(vz0)),
                report,
                ids_name,
                "cylindrical transform of /eqsys/x_p[0] and /eqsys/v_p[0]",
            )

        if not resize_child_aos(injector, "fragment", group_shards.size):
            report.skip(ids_name, f"injector/{igroup}/fragment", "could not resize AoS", "/eqsys/x_p,/eqsys/Y_p")
            continue

        for local_index, ish_raw in enumerate(group_shards):
            ish = int(ish_raw)
            fragment = injector.fragment[local_index]
            report.bind(fragment, "spi.injector.fragment")
            if x_p is not None and ish < x_p.shape[1]:
                r, phi, z = dream_spi_position_to_imas_rphiz(x_p[:, ish, :], R0)
                fill_ids_field(fragment, "position/r", r, report, ids_name, f"cylindrical transform of /eqsys/x_p[:,{ish},:]")
                fill_ids_field(fragment, "position/phi", phi, report, ids_name, f"cylindrical transform of /eqsys/x_p[:,{ish},:]")
                fill_ids_field(fragment, "position/z", z, report, ids_name, f"cylindrical transform of /eqsys/x_p[:,{ish},:]")
            if x_p is not None and v_p is not None and ish < x_p.shape[1] and ish < v_p.shape[1]:
                vr, vphi, vz = dream_spi_velocity_to_imas_rphiz(x_p[:, ish, :], v_p[:, ish, :], R0)
                source = f"cylindrical transform of /eqsys/x_p[:,{ish},:] and /eqsys/v_p[:,{ish},:]"
                fill_ids_field(fragment, "velocity_r", vr, report, ids_name, source)
                fill_ids_field(fragment, "velocity_phi", vphi, report, ids_name, source)
                fill_ids_field(fragment, "velocity_z", vz, report, ids_name, source)
            elif v_p is not None and ish < v_p.shape[1]:
                vr = v_p[:, ish, 0]
                vz = v_p[:, ish, 1]
                fill_ids_field(fragment, "velocity_r", vr, report, ids_name, f"/eqsys/v_p[:,{ish},0]")
                fill_ids_field(fragment, "velocity_phi", np.zeros_like(vr), report, ids_name, "missing /eqsys/x_p; velocity_phi set to zero")
                fill_ids_field(fragment, "velocity_z", vz, report, ids_name, f"/eqsys/v_p[:,{ish},1]")
            if volume is not None and ish < volume.shape[1]:
                fill_ids_field(fragment, "volume", volume[:, ish], report, ids_name, f"/eqsys/Y_p[:,{ish}]")

    return spi


def cell_center_to_edges(xc):
    """
    Construct cell-edge values from cell-center values.

    Interior edges are midpoint averages.
    Boundary edges are linear extrapolations.
    """
    xc = np.asarray(xc)

    xf = np.empty(xc.size + 1, dtype=xc.dtype)

    xf[1:-1] = 0.5 * (xc[:-1] + xc[1:])
    xf[0] = xc[0] - 0.5 * (xc[1] - xc[0])
    xf[-1] = xc[-1] + 0.5 * (xc[-1] - xc[-2])

    return xf


def dVdpsi_from_dVdr(dVdr, dr, psi):
    """
    Convert dV/dr to dV/dpsi while preserving the cell-integrated volume.

    dVdr, dr, psi are all cell-centered arrays with the same length.
    """
    dVdr = np.asarray(dVdr)
    dr = np.asarray(dr)
    psi = np.asarray(psi)

    dV_cell = dVdr * dr

    psi_f = cell_center_to_edges(psi)
    dpsi = np.diff(psi_f)

    dVdpsi = dV_cell / dpsi

    return dVdpsi, dpsi, psi_f


@dataclass
class FluxSurface2DContext:
    psi_on_edges: bool
    theta: np.ndarray
    r_2d: np.ndarray
    z_2d: np.ndarray
    source: str


def prepare_flux_surface_profiles_2d(
    dream: DreamH5,
    R0: float,
    Z0: float,
    psi_size: int,
    report: MappingReport,
) -> Optional[FluxSurface2DContext]:
    Rm = dream.arr("/grid/eq/RMinusR0_f")
    Zm = dream.arr("/grid/eq/ZMinusZ0_f")
    geometry_source = "/grid/eq/RMinusR0_f, /grid/eq/ZMinusZ0_f"
    if Rm is None or Zm is None:
        Rm = dream.arr("/grid/eq/RMinusR0")
        Zm = dream.arr("/grid/eq/ZMinusZ0")
        geometry_source = "/grid/eq/RMinusR0, /grid/eq/ZMinusZ0"
    if Rm is None or Zm is None:
        return None

    Rm = np.asarray(Rm, dtype=float)
    Zm = np.asarray(Zm, dtype=float)
    if Rm.ndim != 2 or Zm.ndim != 2 or Rm.shape != Zm.shape or psi_size == 0:
        report.skip(
            "equilibrium",
            "time_slice/profiles_2d",
            "DREAM flux-surface geometry and psi arrays have incompatible shapes",
            geometry_source,
        )
        return None

    if Rm.shape[1] == psi_size + 1:
        psi_on_edges = True
        psi_source = "cell-edge psi derived from /eqsys/psi_p"
    elif Rm.shape[1] == psi_size:
        psi_on_edges = False
        psi_source = "/eqsys/psi_p"
    else:
        report.skip(
            "equilibrium",
            "time_slice/profiles_2d",
            f"radial geometry size {Rm.shape[1]} does not match psi size {psi_size}",
            f"{geometry_source} and /eqsys/psi_p",
        )
        return None

    theta = flatten_1d(dream.arr("/grid/eq/theta"))
    if theta is None or theta.size != Rm.shape[0]:
        theta = np.linspace(0.0, 2.0 * np.pi, Rm.shape[0], endpoint=False)

    return FluxSurface2DContext(
        psi_on_edges=psi_on_edges,
        theta=theta,
        r_2d=(R0 + Rm).T,
        z_2d=(Z0 + Zm).T,
        source=f"{geometry_source}; {psi_source}",
    )


def flux_surface_profiles_2d(
    dream: DreamH5,
    R0: float,
    Z0: float,
    psi: np.ndarray,
    report: MappingReport,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], str]:
    Rm = dream.arr("/grid/eq/RMinusR0_f")
    Zm = dream.arr("/grid/eq/ZMinusZ0_f")
    geometry_source = "/grid/eq/RMinusR0_f, /grid/eq/ZMinusZ0_f"
    if Rm is None or Zm is None:
        Rm = dream.arr("/grid/eq/RMinusR0")
        Zm = dream.arr("/grid/eq/ZMinusZ0")
        geometry_source = "/grid/eq/RMinusR0, /grid/eq/ZMinusZ0"
    if Rm is None or Zm is None:
        return None, None, None, None, geometry_source

    Rm = np.asarray(Rm, dtype=float)
    Zm = np.asarray(Zm, dtype=float)
    psi = np.asarray(psi, dtype=float).reshape(-1)
    if Rm.ndim != 2 or Zm.ndim != 2 or Rm.shape != Zm.shape or psi.size == 0:
        report.skip(
            "equilibrium",
            "time_slice/profiles_2d",
            "DREAM flux-surface geometry and psi arrays have incompatible shapes",
            geometry_source,
        )
        return None, None, None, None, geometry_source

    if Rm.shape[1] == psi.size + 1:
        psi_grid = cell_center_to_edges(psi)
        psi_source = "cell-edge psi derived from /eqsys/psi_p"
    elif Rm.shape[1] == psi.size:
        psi_grid = psi
        psi_source = "/eqsys/psi_p"
    else:
        report.skip(
            "equilibrium",
            "time_slice/profiles_2d",
            f"radial geometry size {Rm.shape[1]} does not match psi size {psi.size}",
            f"{geometry_source} and /eqsys/psi_p",
        )
        return None, None, None, None, geometry_source

    theta = flatten_1d(dream.arr("/grid/eq/theta"))
    if theta is None or theta.size != Rm.shape[0]:
        theta = np.linspace(0.0, 2.0 * np.pi, Rm.shape[0], endpoint=False)

    r_2d = (R0 + Rm).T
    z_2d = (Z0 + Zm).T

    return psi_grid, theta, r_2d, z_2d, f"{geometry_source}; {psi_source}"


def fill_equilibrium_profiles_2d(
    ts: Any,
    dream: DreamH5,
    R0: float,
    Z0: float,
    psi: np.ndarray,
    report: MappingReport,
    ids_name: str,
    context: Optional[FluxSurface2DContext] = None,
) -> None:
    if not resize_child_aos(ts, "profiles_2d", 1):
        report.skip(ids_name, "time_slice/profiles_2d", "could not resize AoS", "")
        return

    if context is None:
        psi_grid, theta, r_2d, z_2d, source = flux_surface_profiles_2d(dream, R0, Z0, psi, report)
    else:
        psi_grid = cell_center_to_edges(psi) if context.psi_on_edges else np.asarray(psi, dtype=float).reshape(-1)
        theta = context.theta
        r_2d = context.r_2d
        z_2d = context.z_2d
        source = context.source
    if psi_grid is None or theta is None or r_2d is None or z_2d is None:
        return

    p2d = ts.profiles_2d[0]
    ts_path = report.object_paths.get(id(ts), "equilibrium.time_slice")
    report.bind(p2d, f"{ts_path}.profiles_2d")
    fill_ids_field(p2d, "type/name", "total", report, ids_name, "equilibrium_profiles_2d_identifier")
    fill_ids_field(p2d, "type/index", 0, report, ids_name, "equilibrium_profiles_2d_identifier")
    fill_ids_field(p2d, "type/description", "Total fields", report, ids_name, "equilibrium_profiles_2d_identifier")
    fill_ids_field(p2d, "grid_type/name", "inverse_psi_polar", report, ids_name, "poloidal_plane_coordinates_identifier")
    fill_ids_field(p2d, "grid_type/index", 13, report, ids_name, "poloidal_plane_coordinates_identifier")
    fill_ids_field(
        p2d,
        "grid_type/description",
        "Flux-surface grid with psi as radial label and DREAM poloidal angle as dim2",
        report,
        ids_name,
        "poloidal_plane_coordinates_identifier",
    )
    fill_ids_field(p2d, "grid/dim1", psi_grid, report, ids_name, "/eqsys/psi_p")
    fill_ids_field(p2d, "grid/dim2", theta, report, ids_name, "/grid/eq/theta")
    fill_ids_field(p2d, "r", r_2d, report, ids_name, source)
    fill_ids_field(p2d, "z", z_2d, report, ids_name, source)
    fill_ids_field(p2d, "psi", np.repeat(psi_grid[:, np.newaxis], theta.size, axis=1), report, ids_name, source)


def map_equilibrium(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "equilibrium"
    eq = make_ids(factory, ids_name, report)
    if eq is None:
        return None

    time = grids["time"]

    nt = len(time)
    R0 = grids.get("R0") or 0.0
    Z0 = dream.scalar("/grid/eq/Z0", 0.0) or 0.0

    fill_ids_field(eq, "time", time, report, ids_name, "/grid/t")
    fill_vacuum_toroidal_field(eq, ids_name, grids, report)

    if not resize_child_aos(eq, "time_slice", nt):
        report.skip(ids_name, "time_slice", "could not resize AoS", "")
        return eq

    psi_p = scale_optional(time_aligned(dream.arr("/eqsys/psi_p", report), nt), psi_cocos)
    j_tot = time_aligned(dream.arr("/eqsys/j_tot", report), nt) * phi_sign
    ip = flatten_1d(dream.arr("/eqsys/I_p", report)) * phi_sign
    
    r_boundary, z_boundary = boundary_outline(dream, R0, Z0)

    # Some radial grid profiles
    rho_tor = grids["rho_tor"]
    rho_tor_norm = grids["rho_tor_norm"]
    phi_tor = grids["phi_tor"]
    r = grids["r"]
    dr = grids["dr"]

    R_outboard = grids["R_outboard"] 
    dVdr = grids.get("dVdr")
    Bmin = grids.get("Bmin")
    Bmax = grids.get("Bmax")
    FSA_B_over_Bmin = grids.get("FSA_B_over_Bmin")
    FSA_B_over_Bmin2 = grids.get("FSA_B_over_Bmin2")
    R2inv = grids.get("R2inv")
    effective_passing_fraction = grids.get("effective_passing_fraction")

    b_field_average = None
    if Bmin is not None and FSA_B_over_Bmin is not None:
        b_field_average = Bmin * FSA_B_over_Bmin

    gm1 = None
    if R2inv is not None and R0 != 0.0:
        gm1 = R2inv / R0**2

    gm5 = None
    if Bmin is not None and FSA_B_over_Bmin2 is not None:
        gm5 = Bmin**2 * FSA_B_over_Bmin2

    trapped_fraction = None
    trapped_fraction_source = ""
    if effective_passing_fraction is not None:
        trapped_fraction = np.clip(1.0 - effective_passing_fraction, 0.0, 1.0)
        trapped_fraction_source = "derived from /grid/geometry/effectivePassingFraction"

    profiles_2d_context = None
    if psi_p is not None and np.asarray(psi_p).ndim >= 2:
        profiles_2d_context = prepare_flux_surface_profiles_2d(dream, R0, Z0, np.asarray(psi_p).shape[1], report)

    for it in range(nt):
        ts = eq.time_slice[it]
        report.bind(ts, "equilibrium.time_slice")
        fill_ids_field(ts, "time", time[it], report, ids_name, "/grid/t")
        if ip is not None and it < len(ip):
            fill_ids_field(ts, "global_quantities/ip", float(ip[it]), report, ids_name, "/eqsys/I_p")
        
        fill_ids_field(ts, "profiles_1d/rho_tor", rho_tor, report, ids_name, "derived from rho_tor /grid/geometry/toroidalFlux")
        fill_ids_field(ts, "profiles_1d/rho_tor_norm", rho_tor_norm, report, ids_name, "derived from rho_tor")
        fill_ids_field(ts, "profiles_1d/phi", phi_tor, report, ids_name, "derived from phi_tor")
        fill_ids_field(ts, "profiles_1d/r_outboard", R_outboard, report, ids_name, "derived from R0 and r")

        if Bmin is not None:
            fill_ids_field(ts, "profiles_1d/b_field_min", Bmin, report, ids_name, "/grid/geometry/Bmin")
        if Bmax is not None:
            fill_ids_field(ts, "profiles_1d/b_field_max", Bmax, report, ids_name, "/grid/geometry/Bmax")
        if b_field_average is not None:
            fill_ids_field(
                ts,
                "profiles_1d/b_field_average",
                b_field_average,
                report,
                ids_name,
                "derived from /grid/geometry/Bmin and /grid/geometry/FSA_BOverBmin",
            )
        if trapped_fraction is not None:
            fill_ids_field(
                ts,
                "profiles_1d/trapped_fraction",
                trapped_fraction,
                report,
                ids_name,
                trapped_fraction_source,
            )
        if gm1 is not None:
            fill_ids_field(
                ts,
                "profiles_1d/gm1",
                gm1,
                report,
                ids_name,
                "derived from /grid/geometry/FSA_R02OverR2 and /grid/R0",
            )
        if gm5 is not None:
            fill_ids_field(
                ts,
                "profiles_1d/gm5",
                gm5,
                report,
                ids_name,
                "derived from /grid/geometry/Bmin and /grid/geometry/FSA_BOverBmin2",
            )

        if j_tot is not None:
            fill_ids_field(ts, "profiles_1d/j_phi", j_tot[it], report, ids_name, "/eqsys/j_tot")

        if psi_p is not None:
            psi_arr  = np.asarray(psi_p[it], dtype=float)
            psi_bnd  = float(psi_arr[-1]) if psi_arr.size > 0 else 0.0
            psi_axis = float(psi_arr[0]) if psi_arr.size > 0 else 0.0
            psi_norm = normalized_radius(psi_arr) if psi_arr is not None else None

            fill_ids_field(ts, "profiles_1d/psi", psi_arr, report, ids_name, "/eqsys/psi_p")
            fill_ids_field(ts, "profiles_1d/psi_norm", psi_norm, report, ids_name, "derived")

            fill_ids_field(ts, "global_quantities/psi_boundary", psi_bnd, report, ids_name, "/eqsys/psi_p[:,-1]")
            fill_ids_field(ts, "global_quantities/psi_magnetic_axis", psi_axis, report, ids_name, "/eqsys/psi_p[:,0]")

            dVdpsi, dpsi, psi_f = dVdpsi_from_dVdr(dVdr, dr, psi_arr)

            fill_ids_field(ts, "profiles_1d/dvolume_dpsi", dVdpsi, report, ids_name, "derived")
            fill_equilibrium_profiles_2d(ts, dream, R0, Z0, psi_arr, report, ids_name, profiles_2d_context)

        if r_boundary is not None and z_boundary is not None:
            fill_ids_field(ts, "boundary/outline/r", r_boundary, report, ids_name, "/grid/eq/RMinusR0_f[:,-1] + R0")
            fill_ids_field(ts, "boundary/outline/z", z_boundary, report, ids_name, "/grid/eq/ZMinusZ0_f[:,-1] + Z0")
            fill_ids_field(ts, "boundary/geometric_axis/r", float(np.mean([np.min(r_boundary), np.max(r_boundary)])), report, ids_name, "/grid/eq/RMinusR0_f")
            fill_ids_field(ts, "boundary/geometric_axis/z", float(np.mean([np.min(z_boundary), np.max(z_boundary)])), report, ids_name, "/grid/eq/ZMinusZ0_f")

    return eq


def boundary_outline(dream: DreamH5, R0: float, Z0: float) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    Rm = dream.arr("/grid/eq/RMinusR0_f")
    Zm = dream.arr("/grid/eq/ZMinusZ0_f")
    if Rm is None or Zm is None:
        Rm = dream.arr("/grid/eq/RMinusR0")
        Zm = dream.arr("/grid/eq/ZMinusZ0")
    if Rm is None or Zm is None:
        return None, None
    Rm = np.asarray(Rm, dtype=float)
    Zm = np.asarray(Zm, dtype=float)
    if Rm.ndim != 2 or Zm.ndim != 2:
        return None, None
    r_outline = R0 + Rm[:, -1]
    z_outline = Z0 + Zm[:, -1]
    # Close contour if not already closed.
    if r_outline.size > 1 and (r_outline[0] != r_outline[-1] or z_outline[0] != z_outline[-1]):
        r_outline = np.r_[r_outline, r_outline[0]]
        z_outline = np.r_[z_outline, z_outline[0]]
    return r_outline, z_outline


def scalar_time_trace(data: Optional[np.ndarray], nt: int) -> Optional[np.ndarray]:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2 and arr.shape[1] == 1:
        return arr[:, 0]
    return None


def radial_profile_trace(data: Optional[np.ndarray], nt: int) -> Optional[np.ndarray]:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 2:
        return arr
    return None


def radial_sample_trace(data: Optional[np.ndarray], nt: int, index: int) -> Optional[np.ndarray]:
    arr = radial_profile_trace(data, nt)
    if arr is None or arr.shape[1] == 0:
        return None
    return arr[:, index]


def volume_weights(grids: dict[str, Any]) -> Optional[np.ndarray]:
    dVdr = grids.get("dVdr")
    dr = grids.get("dr")
    if dVdr is None or dr is None:
        return None
    weights = np.asarray(dVdr, dtype=float) * np.asarray(dr, dtype=float)
    if weights.ndim != 1 or weights.size == 0 or not np.any(np.isfinite(weights)):
        return None
    return weights


def volume_average_trace(data: Optional[np.ndarray], grids: dict[str, Any], nt: int) -> Optional[np.ndarray]:
    arr = radial_profile_trace(data, nt)
    weights = volume_weights(grids)
    if arr is None or weights is None or arr.shape[1] != weights.size:
        return None
    denom = np.sum(weights)
    if denom == 0:
        return None
    return np.sum(arr * weights[None, :], axis=1) / denom


def volume_integral_trace(data: Optional[np.ndarray], grids: dict[str, Any], nt: int) -> Optional[np.ndarray]:
    aligned = time_aligned(data, nt)
    weights = volume_weights(grids)
    if aligned is None or weights is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == weights.size:
        return np.sum(arr * weights[None, :], axis=1)
    if arr.ndim == 3 and arr.shape[-1] == weights.size:
        return np.sum(arr * weights[None, None, :], axis=(1, 2))
    return None


def trace_time_axis(data: np.ndarray, full_time: np.ndarray) -> Optional[np.ndarray]:
    """Return the time axis matching a DREAM time-dependent output array."""
    arr = np.asarray(data)
    time = np.asarray(full_time, dtype=float)
    if arr.ndim == 0 or time.ndim != 1:
        return None
    if arr.shape[0] == time.size:
        return time
    if arr.shape[0] == time.size - 1:
        return time[1:]
    if arr.shape[0] < time.size:
        return time[: arr.shape[0]]
    return None


def current_from_j_trace(data: Optional[np.ndarray], grids: dict[str, Any], nt: int) -> Optional[np.ndarray]:
    arr = radial_profile_trace(data, nt)
    weight_int_area = grids.get("weight_int_area")
    if arr is None or weight_int_area is None:
        return None
    weights = np.asarray(weight_int_area, dtype=float)
    if weights.ndim != 1 or weights.size != arr.shape[1]:
        return None
    return np.sum(arr * weights[None, :], axis=1)


def fill_code_metadata(ids: Any, dream: DreamH5, report: MappingReport, ids_name: str) -> None:
    """Fill IMAS util/code metadata from DREAM /code provenance datasets."""
    code_values = {
        "commit": dream_text(dream, "/code/commit", report),
        "changes": dream_text(dream, "/code/changes", report),
        "datetime_commit": dream_text(dream, "/code/datetime_commit", report),
        "datetime_simulation": dream_text(dream, "/code/datetime_simulation", report),
        "refspec": dream_text(dream, "/code/refspec", report),
    }
    if not any(code_values.values()):
        return

    fill_ids_field(ids, "code/name", "DREAM", report, ids_name, "DREAM code name")
    fill_ids_field(
        ids,
        "code/description",
        "Disruption Runaway Electron Analysis Model simulation",
        report,
        ids_name,
        "DREAM code description",
    )
    fill_ids_field(
        ids,
        "code/repository",
        "https://github.com/chalmersplasmatheory/DREAM",
        report,
        ids_name,
        "DREAM public repository",
    )

    if code_values["commit"] is not None:
        fill_ids_field(ids, "code/commit", code_values["commit"], report, ids_name, "/code/commit")
    if code_values["refspec"] is not None:
        fill_ids_field(ids, "code/version", code_values["refspec"], report, ids_name, "/code/refspec")

    parameters = {key: value for key, value in code_values.items() if value is not None}
    if parameters:
        fill_ids_field(
            ids,
            "code/parameters",
            dream_code_parameters_xml(parameters),
            report,
            ids_name,
            "/code/changes, /code/commit, /code/datetime_commit, /code/datetime_simulation, /code/refspec",
        )


def map_radiation(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "radiation"
    radiation = make_ids(factory, ids_name, report)
    if radiation is None:
        return None

    emissivity_raw = dream.arr("/other/fluid/Tcold_radiation", report)
    if emissivity_raw is None:
        report.warn(
            "radiation: /other/fluid/Tcold_radiation is absent, so no electron radiation process was filled."
        )
        return radiation

    emissivity = np.asarray(emissivity_raw, dtype=float)
    if emissivity.ndim != 2:
        report.warn(
            f"radiation: /other/fluid/Tcold_radiation has shape {emissivity.shape}, expected (time, radius)."
        )
        return radiation

    time = trace_time_axis(emissivity, grids["time"])
    if time is None:
        report.warn(
            "radiation: could not align /other/fluid/Tcold_radiation to /grid/t, "
            "so no electron radiation process was filled."
        )
        return radiation

    weights = volume_weights(grids)
    if weights is None or weights.size != emissivity.shape[1]:
        report.warn(
            "radiation: volume integral of /other/fluid/Tcold_radiation could not be computed "
            "because /grid/dr, /grid/VpVol, or /grid/R0 do not match the emissivity grid."
        )
        return radiation

    total_power = np.sum(emissivity * weights[None, :], axis=1)
    power_inside = np.cumsum(emissivity * weights[None, :], axis=1)

    fill_ids_field(radiation, "time", time, report, ids_name, "/grid/t aligned to /other/fluid/Tcold_radiation")

    if not resize_child_aos(radiation, "process", 1):
        report.skip(ids_name, "process", "target AoS is not present in this DD version", "/other/fluid/Tcold_radiation")
        return radiation

    process = radiation.process[0]
    report.bind(process, "radiation.process")
    fill_ids_field(process, "identifier/index", -1, report, ids_name, "private process identifier for DREAM total electron radiation")
    fill_ids_field(process, "identifier/name", "dream_total_electron_radiation", report, ids_name, "private process identifier for DREAM total electron radiation")
    fill_ids_field(
        process,
        "identifier/description",
        (
            "Total cold-electron radiated power density from DREAM /other/fluid/Tcold_radiation. "
            "DREAM output provides this as a combined radiation term, so line, recombination, "
            "and bremsstrahlung components are not separated here."
        ),
        report,
        ids_name,
        "/other/fluid/Tcold_radiation",
    )

    if resize_child_aos(process, "global_quantities", time.size):
        for it, t in enumerate(time):
            gq = process.global_quantities[it]
            report.bind(gq, "radiation.process.global_quantities")
            fill_ids_field(gq, "time", float(t), report, ids_name, "/grid/t aligned to /other/fluid/Tcold_radiation")
            fill_ids_field(
                gq,
                "inside_vessel/power_electrons",
                float(total_power[it]),
                report,
                ids_name,
                (
                    "volume integral of /other/fluid/Tcold_radiation using "
                    "dV = /grid/dr * /grid/VpVol * /grid/R0"
                ),
            )
    else:
        report.skip(ids_name, "process/global_quantities", "target AoS is not present in this DD version", "/other/fluid/Tcold_radiation")

    if resize_child_aos(process, "profiles_1d", time.size):
        rho_tor_norm = grids.get("rho_tor_norm")
        rho_tor = grids.get("rho_tor")
        for it, t in enumerate(time):
            p1d = process.profiles_1d[it]
            report.bind(p1d, "radiation.process.profiles_1d")
            fill_ids_field(p1d, "time", float(t), report, ids_name, "/grid/t aligned to /other/fluid/Tcold_radiation")
            if rho_tor_norm is not None:
                fill_ids_field(
                    p1d,
                    "grid/rho_tor_norm",
                    rho_tor_norm,
                    report,
                    ids_name,
                    "rho_tor_norm = normalized sqrt(/grid/geometry/toroidalFlux/(pi*B0))",
                )
            if rho_tor is not None:
                fill_ids_field(
                    p1d,
                    "grid/rho_tor",
                    rho_tor,
                    report,
                    ids_name,
                    "rho_tor = sqrt(/grid/geometry/toroidalFlux/(pi*B0))",
                )
            fill_ids_field(
                p1d,
                "electrons/emissivity",
                emissivity[it],
                report,
                ids_name,
                "/other/fluid/Tcold_radiation radiated power density [J s^-1 m^-3]",
            )
            fill_ids_field(
                p1d,
                "electrons/power_inside",
                power_inside[it],
                report,
                ids_name,
                (
                    "cumulative volume integral of /other/fluid/Tcold_radiation using "
                    "dV = /grid/dr * /grid/VpVol * /grid/R0"
                ),
            )
    else:
        report.skip(ids_name, "process/profiles_1d", "target AoS is not present in this DD version", "/other/fluid/Tcold_radiation")

    return radiation


def map_summary(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "summary"
    summary = make_ids(factory, ids_name, report)
    if summary is None:
        return None
    time = grids["time"]
    nt = len(time)
    fill_ids_field(summary, "time", time, report, ids_name, "/grid/t")
    fill_code_metadata(summary, dream, report, ids_name)

    R0 = grids.get("R0")
    B0 = grids.get("B0")
    if R0 is not None:
        fill_ids_field(summary, "global_quantities/r0/value", float(R0), report, ids_name, "/grid/R0")
        fill_ids_field(summary, "global_quantities/r0/source", "DREAM /grid/R0", report, ids_name, "/grid/R0")
    else:
        report.warn("summary: global_quantities/r0 was not filled because /grid/R0 is missing.")
    if B0 is not None:
        fill_ids_field(summary, "global_quantities/b0/value", np.full(time.shape, float(B0), dtype=float), report, ids_name, "/grid/B0")
        fill_ids_field(summary, "global_quantities/b0/source", "DREAM /grid/B0", report, ids_name, "/grid/B0")
    else:
        report.warn("summary: global_quantities/b0 was not filled because B0 is missing.")

    ip = scalar_time_trace(dream.arr("/eqsys/I_p"), nt) * phi_sign
    n_tot = dream.arr("/eqsys/n_tot")
    n_cold = dream.arr("/eqsys/n_cold")
    n_e_source = "/eqsys/n_tot"
    if n_tot is None:
        n_tot = n_cold
        n_e_source = "/eqsys/n_cold"
    n_re = dream.arr("/eqsys/n_re")
    t_cold = dream.arr("/eqsys/T_cold")
    e_field = dream.arr("/eqsys/E_field")
    zeff = dream.arr("/other/fluid/Zeff")
    j_ohm = dream.arr("/eqsys/j_ohm")  * phi_sign
    j_re = dream.arr("/eqsys/j_re")  * phi_sign
    w_cold = dream.arr("/eqsys/W_cold")
    w_i = dream.arr("/eqsys/W_i")

    if ip is not None:
        fill_ids_field(summary, "global_quantities/ip/value", ip, report, ids_name, "/eqsys/I_p")

    i_ohm = current_from_j_trace(j_ohm, grids, nt)
    i_re = current_from_j_trace(j_re, grids, nt)
    if i_ohm is not None:
        fill_ids_field(summary, "global_quantities/current_ohm/value", i_ohm, report, ids_name, "derived from /eqsys/j_ohm")

    e_cold = volume_integral_trace(w_cold, grids, nt)
    e_ion = volume_integral_trace(w_i, grids, nt)
    e_thermal = None
    if e_cold is not None and e_ion is not None:
        e_thermal = e_cold + e_ion
    elif e_cold is not None:
        e_thermal = e_cold
    if e_cold is not None:
        fill_ids_field(summary, "global_quantities/energy_electrons_thermal/value", e_cold, report, ids_name, "volume integral of /eqsys/W_cold")
    if e_ion is not None:
        fill_ids_field(summary, "global_quantities/energy_ion_total_thermal/value", e_ion, report, ids_name, "volume integral of /eqsys/W_i")
    if e_thermal is not None:
        fill_ids_field(summary, "global_quantities/energy_thermal/value", e_thermal, report, ids_name, "volume integral of /eqsys/W_cold plus /eqsys/W_i when available")

    n_e_volume_average = volume_average_trace(n_tot, grids, nt)
    t_e_volume_average = volume_average_trace(t_cold, grids, nt)
    zeff_volume_average = volume_average_trace(zeff, grids, nt)
    if t_e_volume_average is not None:
        fill_ids_field(summary, "volume_average/t_e/value", t_e_volume_average, report, ids_name, "volume average of /eqsys/T_cold")
    if n_e_volume_average is not None:
        fill_ids_field(summary, "volume_average/n_e/value", n_e_volume_average, report, ids_name, f"volume average of {n_e_source}")
    if zeff_volume_average is not None:
        fill_ids_field(summary, "volume_average/zeff/value", zeff_volume_average, report, ids_name, "volume average of /other/fluid/Zeff")

    local_profiles = {
        "t_e": (t_cold, "/eqsys/T_cold"),
        "n_e": (n_tot, n_e_source),
        "zeff": (zeff, "/other/fluid/Zeff"),
        "e_field_parallel": (e_field, "/eqsys/E_field"),
    }
    for target, (data, source) in local_profiles.items():
        axis_value = radial_sample_trace(data, nt, 0)
        separatrix_value = radial_sample_trace(data, nt, -1)
        if axis_value is not None:
            fill_ids_field(summary, f"local/magnetic_axis/{target}/value", axis_value, report, ids_name, f"{source}[:,0]")
        if separatrix_value is not None:
            fill_ids_field(summary, f"local/separatrix/{target}/value", separatrix_value, report, ids_name, f"{source}[:,-1]")


    if i_re is not None:
        fill_ids_field(summary, "runaways/current/value", i_re, report, ids_name, "derived from /eqsys/j_re")
    if i_re is not None and i_re.size > 0:
        fill_ids_field(summary, "runaways/current_phi_max/value", float(np.nanmax(np.abs(i_re))), report, ids_name, "maximum absolute derived runaway current")
    runaway_particles = volume_integral_trace(n_re, grids, nt)
    if runaway_particles is not None:
        fill_ids_field(summary, "runaways/particles/value", runaway_particles, report, ids_name, "volume integral of /eqsys/n_re")

    return summary


# ----------------------------- writing ------------------------------------

def write_ids(ids_list: list[Any], uri: str, report: MappingReport) -> None:
    ensure_imas()
    db = imas.DBEntry(uri, "w")
    try:
        for ids in ids_list:
            if ids is not None:
                db.put(ids)
        report.written_uri = uri
    finally:
        try:
            db.close()
        except Exception:
            pass


def timed_map(report: MappingReport, key: str, mapper: Any, *args: Any) -> Any:
    start = time.perf_counter()
    ids = mapper(*args)
    report.add_timing(key, time.perf_counter() - start)
    return ids


def build_ids(dream_file: str, dd_version: str | None, selected: Iterable[str]) -> tuple[list[Any], MappingReport]:
    factory = make_factory(dd_version)
    report = MappingReport(source_file=dream_file, dd_version_requested=dd_version)
    dream = DreamH5(dream_file)
    try:
        grids = common_grids(dream, report)
        ids_list = []
        selected_set = set(selected)
        if "plasma_profiles" in selected_set:
            ids_list.append(timed_map(report, "plasma_profiles", map_plasma_profiles, factory, dream, grids, report))
        if "runaway_electrons" in selected_set:
            ids_list.append(timed_map(report, "runaway_electrons", map_runaway_electrons, factory, dream, grids, report))
        if "spi" in selected_set:
            ids_list.append(timed_map(report, "spi", map_spi, factory, dream, grids, report))
        if "equilibrium" in selected_set:
            ids_list.append(timed_map(report, "equilibrium", map_equilibrium, factory, dream, grids, report))
        if "radiation" in selected_set:
            ids_list.append(timed_map(report, "radiation", map_radiation, factory, dream, grids, report))
        if "summary" in selected_set:
            ids_list.append(timed_map(report, "summary", map_summary, factory, dream, grids, report))
        return [ids for ids in ids_list if ids is not None], report
    finally:
        dream.close()


def convert_dream_h5(
    dream_file: str | Path,
    dd_version: str | None = None,
    ids: Iterable[str] | None = None,
) -> ConversionResult:
    """Convert one DREAM HDF5 file into IMAS IDS objects and a mapping report."""
    selected = DEFAULT_IDS if ids is None else list(ids)
    ids_list, report = build_ids(str(dream_file), dd_version, selected)
    return ConversionResult(source_file=str(dream_file), ids_list=ids_list, report=report)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert selected DREAM HDF5 output quantities to IMAS IDSs."
    )
    parser.add_argument("dream_h5", help="Input DREAM HDF5 output file")
    parser.add_argument(
        "--uri",
        default="dream_imas.nc",
        help=(
            "Output URI passed to imas.DBEntry. Examples: 'dream_imas.nc' for netCDF, "
            "or 'imas:hdf5?path=./imas_dream' for IMAS-Core HDF5."
        ),
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        default=DEFAULT_IDS,
        choices=DEFAULT_IDS,
        help="IDSs to create/write.",
    )
    parser.add_argument(
        "--dd-version",
        default=None,
        help="Optional IMAS Data Dictionary version, e.g. 4.1.0. Default: environment/latest.",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Path to write mapping report. Default: <uri or input stem>.mapping_report.md",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Build IDS objects and write only the mapping report; do not call DBEntry.put().",
    )
    return parser.parse_args(argv)


def default_report_path(args: argparse.Namespace) -> Path:
    if args.report:
        return Path(args.report)
    uri_name = args.uri.replace(":", "_").replace("?", "_").replace("/", "_").replace(";", "_")
    if uri_name:
        return Path(f"{uri_name}.mapping_report.md")
    return Path(args.dream_h5).with_suffix(".mapping_report.md")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if not Path(args.dream_h5).exists():
        print(f"ERROR: input file does not exist: {args.dream_h5}", file=sys.stderr)
        return 2
    ids_list, report = build_ids(args.dream_h5, args.dd_version, args.ids)
    if not args.dry_run:
        write_ids(ids_list, args.uri, report)
    report_path = default_report_path(args)
    report.write(report_path)
    print(f"Built IDSs: {', '.join(report.created_ids) if report.created_ids else '(none)'}")
    if args.dry_run:
        print("Dry run: no IMAS DBEntry was written.")
    else:
        print(f"Wrote IMAS data entry: {args.uri}")
    print(f"Mapping report: {report_path}")
    print(
        "Set mappings: "
        f"{len(report.set_nodes)} unique / {sum(report.set_nodes.values())} calls | "
        f"Skipped: {len(report.skipped_nodes)} unique / {sum(report.skipped_nodes.values())} calls | "
        f"Missing sources: {len(report.missing_sources)} unique / {sum(report.missing_sources.values())} calls"
    )
    if report.warnings:
        print("Warnings:")
        warning_items = sorted(report.warnings.items(), key=lambda item: (-item[1], item[0]))
        for w, count in warning_items[:10]:
            suffix = f" ({count}x)" if count > 1 else ""
            print(f"  - {w}{suffix}")
        if len(report.warnings) > 10:
            print(f"  ... {len(report.warnings) - 10} more warnings in the report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
