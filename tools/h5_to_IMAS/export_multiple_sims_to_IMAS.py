#!/usr/bin/env python3
"""
Export multi-stage DREAM simulations to one IMAS entry per simulation.

Each DREAM stage is converted with dh5file_to_imas.py, shifted to a global
simulation time, overlap-trimmed, and written with DBEntry.put()/put_slice().
"""

from __future__ import annotations

import argparse
import copy
import getpass
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import h5py
import numpy as np

from h5file_to_imas import (
    DEFAULT_IDS,
    ConversionResult,
    MappingReport,
    convert_dream_h5,
    ensure_imas,
    imas,
    resize_aos,
    fill_ids_field,
)


TIME_TOLERANCE = 1e-12
DEFAULT_URI_TEMPLATE = "imas:hdf5?user={user};pulse={pulse};run={run};database={database};version={version}"
STAGE_ORDER = {
    "init": 0,
    "injection_start": 10,
    "injection_continuation": 20,
    "injection": 25,
    "cq": 30,
    "unknown": 100,
}
COMMON_TOP_LEVEL_TIME_FIELDS = ("vacuum_toroidal_field/b0",)
TOP_LEVEL_TIME_FIELDS = {
    "plasma_profiles": COMMON_TOP_LEVEL_TIME_FIELDS,
    "runaway_electrons": COMMON_TOP_LEVEL_TIME_FIELDS + ("global_quantities/current_phi",),
    "equilibrium": COMMON_TOP_LEVEL_TIME_FIELDS,
    "summary": (
        "global_quantities/ip/value",
        "global_quantities/current_ohm/value",
        "global_quantities/energy_electrons_thermal/value",
        "global_quantities/energy_ion_total_thermal/value",
        "global_quantities/energy_thermal/value",
        "volume_average/t_e/value",
        "volume_average/n_e/value",
        "volume_average/zeff/value",
        "local/magnetic_axis/t_e/value",
        "local/magnetic_axis/n_e/value",
        "local/magnetic_axis/zeff/value",
        "local/magnetic_axis/e_field_parallel/value",
        "local/separatrix/t_e/value",
        "local/separatrix/n_e/value",
        "local/separatrix/zeff/value",
        "local/separatrix/e_field_parallel/value",
        "runaways/current/value",
        "runaways/particles/value",
    ),
}


@dataclass
class StageResult:
    path: Path
    kind: str
    conversion: ConversionResult
    local_time: np.ndarray
    global_time: np.ndarray
    offset: float
    previous_end: Optional[float]
    read_seconds: float = 0.0
    convert_seconds: float = 0.0


@dataclass
class SimulationResult:
    name: str
    directory: Path
    uri: str
    pulse: int
    run: int
    stages: list[StageResult] = field(default_factory=list)
    report: Optional[MappingReport] = None
    timings: dict[str, float] = field(default_factory=dict)


def natural_text_key(text: str) -> list[Any]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def natural_key(path: Path) -> list[Any]:
    return natural_text_key(path.name)


def natural_relative_key(path: Path, root: Path) -> list[Any]:
    try:
        text = path.relative_to(root).as_posix()
    except ValueError:
        text = path.as_posix()
    return natural_text_key(text)


def classify_stage(path: Path) -> str:
    name = path.stem.lower()
    if re.search(r"(^|_)init(ial)?($|_)", name):
        return "init"
    if "injection_start" in name:
        return "injection_start"
    if "injection_continuation" in name:
        return "injection_continuation"
    if "injection" in name:
        return "injection"
    if re.search(r"(^|_)cq($|_)", name) or "current_quench" in name:
        return "cq"
    return "unknown"


def sort_stage_files(files: Iterable[Path]) -> list[Path]:
    return sorted(files, key=lambda path: (STAGE_ORDER[classify_stage(path)], natural_key(path)))


def keep_stage_file(path: Path) -> bool:
    return classify_stage(path) != "init"


def discover_simulations(root: Path, stage_glob: str) -> list[tuple[str, Path, list[Path]]]:
    simulations: list[tuple[str, Path, list[Path]]] = []

    output_dirs = {
        path.parent
        for path in root.rglob(stage_glob)
        if path.is_file() and path.parent.name == "output"
    }
    for output_dir in sorted(output_dirs, key=lambda path: natural_relative_key(path.parent, root)):
        directory = output_dir.parent
        files = [path for path in sort_stage_files(output_dir.glob(stage_glob)) if path.is_file() and keep_stage_file(path)]
        if files:
            simulations.append((simulation_name(root, directory), directory, files))
    if simulations:
        return simulations

    direct_files = [path for path in sort_stage_files(root.glob(stage_glob)) if path.is_file() and keep_stage_file(path)]
    if direct_files:
        return [(root.name, root, direct_files)]

    for child in sorted((path for path in root.iterdir() if path.is_dir()), key=lambda path: natural_relative_key(path, root)):
        files = [path for path in sort_stage_files(child.glob(stage_glob)) if path.is_file() and keep_stage_file(path)]
        if files:
            simulations.append((simulation_name(root, child), child, files))
    return simulations


def simulation_name(root: Path, directory: Path) -> str:
    try:
        return directory.relative_to(root).as_posix()
    except ValueError:
        return directory.name


def read_dream_time(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h5:
        if "/grid/t" not in h5:
            raise ValueError(f"{path} does not contain /grid/t")
        return np.asarray(h5["/grid/t"][()], dtype=float).reshape(-1)


def stage_global_time(local_time: np.ndarray, previous_end: Optional[float]) -> tuple[np.ndarray, float]:
    if local_time.size == 0 or previous_end is None:
        return local_time, 0.0
    offset = previous_end - float(local_time[0]) if local_time[0] <= previous_end + TIME_TOLERANCE else 0.0
    return local_time + offset, offset


def keep_after_overlap(local_time: np.ndarray, offset: float, previous_end: Optional[float]) -> tuple[np.ndarray, np.ndarray]:
    global_time = np.asarray(local_time, dtype=float) + offset
    keep = np.ones(global_time.shape, dtype=bool) if previous_end is None else global_time > previous_end + TIME_TOLERANCE
    return global_time, keep


def ids_by_name(conversion: ConversionResult) -> dict[str, Any]:
    return dict(zip(conversion.report.created_ids, conversion.ids_list))


def merge_reports(source_name: str, reports: Iterable[MappingReport]) -> MappingReport:
    merged = MappingReport(source_file=source_name)
    for report in reports:
        merged.dd_version_requested = merged.dd_version_requested or report.dd_version_requested
        for ids_name in report.created_ids:
            if ids_name not in merged.created_ids:
                merged.created_ids.append(ids_name)
        merge_counts(merged.set_nodes, report.set_nodes)
        merge_counts(merged.skipped_nodes, report.skipped_nodes)
        merge_counts(merged.missing_sources, report.missing_sources)
        merge_counts(merged.warnings, report.warnings)
        for key, seconds in getattr(report, "timings", {}).items():
            merged.add_timing(key, seconds)
    return merged


def merge_counts(target: dict[Any, int], source: dict[Any, int]) -> None:
    for key, count in source.items():
        target[key] = target.get(key, 0) + count


def get_time(node: Any) -> Optional[np.ndarray]:
    if not hasattr(node, "time"):
        return None
    try:
        return np.asarray(getattr(node, "time"), dtype=float)
    except Exception:
        return None


def get_node(root: Any, path: str) -> Any:
    node = root
    for part in [item for item in path.split("/") if item]:
        node = node[int(part)] if part.isdigit() else getattr(node, part)
    return node


def set_stage_time(node: Any, time: np.ndarray, ids_name: str, report: MappingReport) -> None:
    fill_ids_field(node, "time", np.asarray(time, dtype=float), report, ids_name, "stage time offset and overlap trim for put_slice")


def trim_numeric_field(root: Any, path: str, keep: np.ndarray, original_length: int, ids_name: str, report: MappingReport) -> None:
    parts = [part for part in path.split("/") if part]
    if not parts:
        return
    try:
        parent = get_node(root, "/".join(parts[:-1])) if len(parts) > 1 else root
        value = getattr(parent, parts[-1])
        arr = np.asarray(value, dtype=float)
    except Exception:
        return
    if arr.ndim < 1 or arr.shape[0] != original_length:
        return
    if len(parts) > 1:
        report.bind(parent, f"{ids_name}.{'.'.join(parts[:-1])}")
    fill_ids_field(
        parent,
        parts[-1],
        arr[keep],
        report,
        ids_name,
        f"{path} trimmed to stage time for put_slice",
    )


def trim_top_level_time_fields(ids_name: str, ids: Any, keep: np.ndarray, original_length: int, report: MappingReport) -> None:
    for path in TOP_LEVEL_TIME_FIELDS.get(ids_name, ()):
        trim_numeric_field(ids, path, keep, original_length, ids_name, report)


def trim_numeric_traces(node: Any, keep: np.ndarray, original_length: int, seen: Optional[set[int]] = None) -> None:
    if seen is None:
        seen = set()
    if id(node) in seen:
        return
    seen.add(id(node))

    for name in dir(node):
        if name.startswith("_") or name in {"metadata", "coordinates", "time", "copy_from", "validate"}:
            continue
        try:
            value = getattr(node, name)
        except Exception:
            continue
        if callable(value) or (hasattr(value, "resize") and hasattr(value, "__getitem__")):
            continue
        try:
            arr = np.asarray(value, dtype=float)
        except Exception:
            try:
                trim_numeric_traces(value, keep, original_length, seen)
            except Exception:
                pass
            continue
        if arr.ndim >= 1 and arr.shape[0] == original_length:
            try:
                setattr(node, name, arr[keep])
            except Exception:
                pass


def trim_aos(parent: Any, aos_name: str, keep: np.ndarray, item_times: Optional[np.ndarray], ids_name: str, report: MappingReport) -> None:
    if not hasattr(parent, aos_name):
        return
    aos = getattr(parent, aos_name)
    indices = np.flatnonzero(keep)
    if indices.size == len(aos) and np.array_equal(indices, np.arange(len(aos))):
        if item_times is not None:
            for index in indices:
                if hasattr(aos[int(index)], "time"):
                    fill_ids_field(aos[int(index)], "time", float(item_times[int(index)]), report, ids_name, "stage time offset for put_slice")
        return

    for out_index, source_index in enumerate(indices):
        if out_index == int(source_index):
            continue
        try:
            aos[out_index].copy_from(aos[int(source_index)])
        except Exception:
            try:
                aos[out_index] = copy.deepcopy(aos[int(source_index)])
            except Exception as exc:
                report.skip(ids_name, aos_name, f"could not copy trimmed AoS item: {type(exc).__name__}: {exc}", "stage trim for put_slice")

    if not resize_aos(aos, len(indices)):
        report.skip(ids_name, aos_name, "could not trim target AoS", "stage trim for put_slice")
        return
    for out_index in range(len(indices)):
        if item_times is not None and hasattr(aos[out_index], "time"):
            fill_ids_field(aos[out_index], "time", float(item_times[indices[out_index]]), report, ids_name, "stage time offset for put_slice")


def trim_standard_ids(
    ids_name: str,
    ids: Any,
    aos_name: Optional[str],
    offset: float,
    previous_end: Optional[float],
    report: MappingReport,
) -> bool:
    time = get_time(ids)
    if time is None:
        report.warn(f"{ids_name}: no top-level time vector found; stage skipped.")
        return False
    global_time, keep = keep_after_overlap(time, offset, previous_end)
    if not np.any(keep):
        report.warn(f"{ids_name}: stage has no time samples after overlap trimming; skipped.")
        return False
    set_stage_time(ids, global_time[keep], ids_name, report)
    trim_top_level_time_fields(ids_name, ids, keep, len(time), report)
    if aos_name is not None:
        trim_aos(ids, aos_name, keep, global_time, ids_name, report)
    return True


def trim_radiation_ids(ids: Any, offset: float, previous_end: Optional[float], report: MappingReport) -> bool:
    ids_name = "radiation"
    wrote_any = False
    time = get_time(ids)
    if time is not None:
        global_time, keep = keep_after_overlap(time, offset, previous_end)
        if np.any(keep):
            set_stage_time(ids, global_time[keep], ids_name, report)
            wrote_any = True

    if not hasattr(ids, "process"):
        return wrote_any
    for iprocess in range(len(ids.process)):
        process = ids.process[iprocess]
        for aos_name in ("global_quantities", "profiles_1d"):
            if not hasattr(process, aos_name):
                continue
            aos = getattr(process, aos_name)
            item_times = []
            for i in range(len(aos)):
                if not hasattr(aos[i], "time"):
                    item_times = []
                    break
                item_times.append(float(np.asarray(getattr(aos[i], "time"), dtype=float)))
            if not item_times:
                continue
            global_times, keep = keep_after_overlap(np.asarray(item_times), offset, previous_end)
            trim_aos(process, aos_name, keep, global_times, ids_name, report)
            wrote_any = wrote_any or np.any(keep)
    return wrote_any


def trim_spi_ids(ids: Any, offset: float, previous_end: Optional[float], report: MappingReport) -> bool:
    ids_name = "spi"
    if not hasattr(ids, "injector") or len(ids.injector) == 0:
        report.warn("spi: stage has no injector content; skipped.")
        return False

    time = get_time(ids)
    if time is None:
        report.warn("spi: no top-level time vector found; stage skipped.")
        return False
    global_time, keep = keep_after_overlap(time, offset, previous_end)
    if not np.any(keep):
        report.warn("spi: stage has no time samples after overlap trimming; skipped.")
        return False

    set_stage_time(ids, global_time[keep], ids_name, report)
    original_length = len(time)

    for injector_index in range(len(ids.injector)):
        injector = ids.injector[injector_index]
        if not hasattr(injector, "fragment"):
            continue
        for fragment_index in range(len(injector.fragment)):
            fragment = injector.fragment[fragment_index]
            report.bind(fragment, f"spi.injector.fragment")
            for path in (
                "position/r",
                "position/z",
                "position/phi",
                "velocity_r",
                "velocity_z",
                "velocity_phi",
                "volume",
            ):
                trim_numeric_field(fragment, path, keep, original_length, ids_name, report)

    return True


def numeric_field(root: Any, path: str) -> Optional[np.ndarray]:
    try:
        value = get_node(root, path)
        arr = np.asarray(value, dtype=float)
    except Exception:
        return None
    return arr if arr.size > 0 else None


def append_numeric_field(root: Any, path: str, values: np.ndarray, report: MappingReport, ids_name: str, source: str) -> None:
    current = numeric_field(root, path)
    if current is None:
        return
    fill_ids_field(root, path, np.concatenate([current, values]), report, ids_name, source)


def merge_spi_stage(base: Any, stage: Any, report: MappingReport) -> None:
    base_time = get_time(base)
    stage_time = get_time(stage)
    if base_time is None or stage_time is None:
        return
    fill_ids_field(base, "time", np.concatenate([base_time, stage_time]), report, "spi", "SPI stages concatenated in memory")

    n_injectors = min(len(base.injector), len(stage.injector)) if hasattr(base, "injector") and hasattr(stage, "injector") else 0
    for injector_index in range(n_injectors):
        base_injector = base.injector[injector_index]
        stage_injector = stage.injector[injector_index]
        if not hasattr(base_injector, "fragment") or not hasattr(stage_injector, "fragment"):
            continue
        n_fragments = min(len(base_injector.fragment), len(stage_injector.fragment))
        for fragment_index in range(n_fragments):
            base_fragment = base_injector.fragment[fragment_index]
            stage_fragment = stage_injector.fragment[fragment_index]
            report.bind(base_fragment, "spi.injector.fragment")
            for path in (
                "position/r",
                "position/z",
                "position/phi",
                "velocity_r",
                "velocity_z",
                "velocity_phi",
                "volume",
            ):
                values = numeric_field(stage_fragment, path)
                if values is not None:
                    append_numeric_field(base_fragment, path, values, report, "spi", "SPI stages concatenated in memory")


def merged_spi_for_write(stages: list[StageResult], report: MappingReport) -> Optional[Any]:
    merged = None
    for stage in stages:
        spi = ids_by_name(stage.conversion).get("spi")
        if spi is None:
            continue
        if not trim_spi_ids(spi, stage.offset, stage.previous_end, report):
            continue
        if merged is None:
            merged = spi
        else:
            merge_spi_stage(merged, spi, report)
    return merged


def prepare_ids_for_stage_write(ids_name: str, ids: Any, offset: float, previous_end: Optional[float], report: MappingReport) -> bool:
    if ids_name in {"plasma_profiles", "runaway_electrons"}:
        return trim_standard_ids(ids_name, ids, "profiles_1d", offset, previous_end, report)
    if ids_name == "equilibrium":
        return trim_standard_ids(ids_name, ids, "time_slice", offset, previous_end, report)
    if ids_name == "summary":
        return trim_standard_ids(ids_name, ids, None, offset, previous_end, report)
    if ids_name == "radiation":
        return trim_radiation_ids(ids, offset, previous_end, report)
    if ids_name == "spi":
        return trim_spi_ids(ids, offset, previous_end, report)

    time = get_time(ids)
    if time is None:
        return previous_end is None
    global_time, keep = keep_after_overlap(time, offset, previous_end)
    if not np.any(keep):
        return False
    set_stage_time(ids, global_time[keep], ids_name, report)
    trim_numeric_traces(ids, keep, len(time))
    return True


def propagate_global_summary_scalars(stages: list[StageResult], report: MappingReport) -> None:
    path = "runaways/current_phi_max/value"
    values = []
    for stage in stages:
        summary = ids_by_name(stage.conversion).get("summary")
        if summary is None:
            continue
        try:
            values.append(float(np.asarray(get_node(summary, path), dtype=float).reshape(-1)[0]))
        except Exception:
            pass
    if not values:
        return
    for stage in stages:
        summary = ids_by_name(stage.conversion).get("summary")
        if summary is not None:
            fill_ids_field(summary, path, max(values), report, "summary", "maximum over all DREAM stages before DBEntry.put()/put_slice()")


def add_timing(timings: dict[str, float], key: str, seconds: float) -> None:
    timings[key] = timings.get(key, 0.0) + seconds


def write_stages_with_put_slice(
    stages: list[StageResult],
    uri: str,
    selected_ids: Iterable[str],
    report: MappingReport,
    profile: bool = False,
) -> dict[str, float]:
    ensure_imas()
    db = imas.DBEntry(uri, "w")
    written_ids: set[str] = set()
    timings: dict[str, float] = {}
    try:
        start = time.perf_counter()
        propagate_global_summary_scalars(stages, report)
        add_timing(timings, "write.propagate_summary", time.perf_counter() - start)
        selected_ids = list(selected_ids)
        if "spi" in selected_ids:
            start = time.perf_counter()
            spi = merged_spi_for_write(stages, report)
            prepare_seconds = time.perf_counter() - start
            add_timing(timings, "write.prepare.spi", prepare_seconds)
            if spi is not None:
                start = time.perf_counter()
                db.put(spi)
                written_ids.add("spi")
                put_seconds = time.perf_counter() - start
                add_timing(timings, "write.put.spi", put_seconds)
                report.warn("spi: concatenated DREAM stages in memory and wrote once with DBEntry.put().")
                if profile:
                    print(f"    spi: prepare/merge {prepare_seconds:.3f}s, put {put_seconds:.3f}s")

        for stage in stages:
            stage_ids = ids_by_name(stage.conversion)
            for ids_name in selected_ids:
                if ids_name == "spi":
                    continue
                ids = stage_ids.get(ids_name)
                if ids is None:
                    continue
                start = time.perf_counter()
                if not prepare_ids_for_stage_write(ids_name, ids, stage.offset, stage.previous_end, report):
                    add_timing(timings, f"write.prepare.{ids_name}", time.perf_counter() - start)
                    continue
                prepare_seconds = time.perf_counter() - start
                add_timing(timings, f"write.prepare.{ids_name}", prepare_seconds)
                start = time.perf_counter()
                if ids_name in written_ids:
                    db.put_slice(ids)
                    write_action = "put_slice"
                    report.warn(f"{ids_name}: appended stage {stage.path.name} with DBEntry.put_slice().")
                else:
                    db.put(ids)
                    written_ids.add(ids_name)
                    write_action = "put"
                    report.warn(f"{ids_name}: wrote first available stage {stage.path.name} with DBEntry.put().")
                write_seconds = time.perf_counter() - start
                add_timing(timings, f"write.{write_action}.{ids_name}", write_seconds)
                if profile:
                    print(
                        f"    {ids_name}: prepare {prepare_seconds:.3f}s, "
                        f"{write_action} {write_seconds:.3f}s"
                    )
        report.written_uri = uri
        return timings
    finally:
        try:
            db.close()
        except Exception:
            pass


def simulation_uri(name: str, directory: Path, args: argparse.Namespace, index: int) -> tuple[str, int, int]:
    pulse = args.pulse_start + index - 1
    if pulse > args.pulse_end:
        raise ValueError(
            f"pulse {pulse} exceeds --pulse-end {args.pulse_end}; "
            "increase the pulse range or export fewer simulations."
        )
    run = args.run
    uri = args.uri_template.format(
        sim=name,
        case=name,
        sim_dir=str(directory),
        index=index,
        pulse=pulse,
        run=run,
        user=args.imas_user,
        database=args.database,
        version=args.version,
    )
    return uri, pulse, run


def print_profile_summary(result: SimulationResult) -> None:
    if not result.timings:
        return
    print("  profile:")
    for key, seconds in sorted(result.timings.items(), key=lambda item: item[1], reverse=True):
        print(f"    {key}: {seconds:.3f}s")


def convert_simulation(name: str, directory: Path, files: list[Path], args: argparse.Namespace, index: int) -> SimulationResult:
    uri, pulse, run = simulation_uri(name, directory, args, index)
    result = SimulationResult(name=name, directory=directory, uri=uri, pulse=pulse, run=run)
    previous_end: Optional[float] = None

    for file_path in files:
        start = time.perf_counter()
        local_time = read_dream_time(file_path)
        read_seconds = time.perf_counter() - start
        global_time, offset = stage_global_time(local_time, previous_end)
        start = time.perf_counter()
        conversion = convert_dream_h5(file_path, args.dd_version, args.ids)
        convert_seconds = time.perf_counter() - start
        result.stages.append(StageResult(
            path=file_path,
            kind=classify_stage(file_path),
            conversion=conversion,
            local_time=local_time,
            global_time=global_time,
            offset=offset,
            previous_end=previous_end,
            read_seconds=read_seconds,
            convert_seconds=convert_seconds,
        ))
        add_timing(result.timings, "stage.read_time", read_seconds)
        add_timing(result.timings, "stage.convert", convert_seconds)
        if args.profile:
            print(f"  {file_path.name}: read time {read_seconds:.3f}s, convert {convert_seconds:.3f}s")
        if global_time.size > 0:
            previous_end = float(global_time[-1])

    result.report = merge_reports(name, (stage.conversion.report for stage in result.stages))
    for key, seconds in result.report.timings.items():
        add_timing(result.timings, f"convert.{key}", seconds)
    result.report.written_uri = None if args.dry_run else uri
    for stage in result.stages:
        if stage.kind == "unknown":
            result.report.warn(f"{name}: stage file {stage.path.name} did not match known DREAM stage names.")
    if not args.dry_run:
        write_timings = write_stages_with_put_slice(result.stages, uri, args.ids, result.report, args.profile)
        for key, seconds in write_timings.items():
            add_timing(result.timings, key, seconds)
    if args.profile:
        print_profile_summary(result)
    return result


def time_range_text(time: np.ndarray) -> str:
    if time.size == 0:
        return "(empty)"
    return f"{float(time[0]):.9g} to {float(time[-1]):.9g} s ({time.size})"


def write_batch_report(results: list[SimulationResult], path: Path) -> None:
    lines = ["# DREAM Multi-Stage IMAS Export Report", ""]
    for result in results:
        lines.extend([
            f"## {result.name}",
            "",
            f"- Directory: `{result.directory}`",
            f"- URI: `{result.uri}`",
            f"- Pulse/run: `{result.pulse}/{result.run}`",
            f"- IDSs: {', '.join(f'`{name}`' for name in result.report.created_ids) if result.report else '(none)'}",
            "",
            "| Stage | Type | File | Local Time | Offset | Global Time |",
            "| --- | --- | --- | --- | --- | --- |",
        ])
        for i, stage in enumerate(result.stages, start=1):
            lines.append(
                f"| {i} | `{stage.kind}` | `{stage.path.name}` | "
                f"{time_range_text(stage.local_time)} | {stage.offset:.9g} | {time_range_text(stage.global_time)} |"
            )
        lines.append("")
        if result.report is not None:
            lines.extend([
                f"- Set mappings: {len(result.report.set_nodes)} unique, {sum(result.report.set_nodes.values())} total calls",
                f"- Skipped mappings: {len(result.report.skipped_nodes)} unique, {sum(result.report.skipped_nodes.values())} total calls",
                f"- Missing sources: {len(result.report.missing_sources)} unique, {sum(result.report.missing_sources.values())} total calls",
                f"- Warnings: {len(result.report.warnings)} unique, {sum(result.report.warnings.values())} total calls",
                "",
            ])
        if result.timings:
            lines.extend([
                "| Timing Section | Seconds |",
                "| --- | --- |",
            ])
            for key, seconds in sorted(result.timings.items(), key=lambda item: item[1], reverse=True):
                lines.append(f"| `{key}` | {seconds:.3f} |")
            lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_uri_map(results: list[SimulationResult], path: Path) -> None:
    lines = ["URI\tProvenance folder"]
    lines.extend(f"{result.uri}\t{result.name}" for result in results)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export multi-stage DREAM simulations to IMAS.")
    parser.add_argument("root", help="Simulation folder, or a folder containing simulation subfolders.")
    parser.add_argument("--stage-glob", default="output*.h5", help="Glob used to find DREAM output stages.")
    parser.add_argument(
        "--uri-template",
        default=DEFAULT_URI_TEMPLATE,
        help=(
            "IMAS URI template. Available fields: {sim}, {case}, {sim_dir}, {index}, "
            "{pulse}, {run}, {user}, {database}, {version}."
        ),
    )
    parser.add_argument("--imas-user", default=getpass.getuser(), help="IMAS user field used by the default URI template.")
    parser.add_argument("--database", default="ITER_DISRUPTIONS_DREAM", help="IMAS database field used by the default URI template.")
    parser.add_argument("--version", default="4", help="IMAS version field used by the default URI template.")
    parser.add_argument("--pulse-start", type=int, default=120000, help="First pulse number assigned to discovered simulations.")
    parser.add_argument("--pulse-end", type=int, default=129999, help="Maximum pulse number allowed for this export.")
    parser.add_argument("--run", type=int, default=1, help="Run number assigned to each exported simulation.")
    parser.add_argument("--ids", nargs="+", default=DEFAULT_IDS, choices=DEFAULT_IDS, help="IDSs to create/write.")
    parser.add_argument("--dd-version", default=None, help="Optional IMAS Data Dictionary version.")
    parser.add_argument("--report", default="dream_multi_stage_export_report.md", help="Markdown batch report path.")
    parser.add_argument("--uri-map", default="dream_simulation_uris.tsv", help="TSV file listing URI and provenance folder.")
    parser.add_argument("--write-simulation-reports", action="store_true", help="Also write one report per simulation.")
    parser.add_argument("--dry-run", action="store_true", help="Convert stages and write reports, but do not write IMAS entries.")
    parser.add_argument("--profile", action="store_true", help="Print and report coarse timing for conversion and IMAS writes.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        print(f"ERROR: root path does not exist: {root}", file=sys.stderr)
        return 2

    simulations = discover_simulations(root, args.stage_glob)
    if not simulations:
        print(f"ERROR: no stage files matching {args.stage_glob!r} found under {root}", file=sys.stderr)
        return 2
    if args.pulse_start > args.pulse_end:
        print(f"ERROR: --pulse-start {args.pulse_start} is larger than --pulse-end {args.pulse_end}", file=sys.stderr)
        return 2
    if args.pulse_start + len(simulations) - 1 > args.pulse_end:
        print(
            f"ERROR: {len(simulations)} simulations need pulses through {args.pulse_start + len(simulations) - 1}, "
            f"but --pulse-end is {args.pulse_end}",
            file=sys.stderr,
        )
        return 2

    results = []
    for index, (name, directory, files) in enumerate(simulations, start=1):
        print(f"[{index}/{len(simulations)}] {name}: {len(files)} stage file(s)")
        print("  stage order: " + " -> ".join(f"{classify_stage(path)}:{path.name}" for path in files))
        result = convert_simulation(name, directory, files, args, index)
        results.append(result)
        if args.write_simulation_reports and result.report is not None:
            report_path = directory / "dream_to_imas_mapping_report.md"
            result.report.write(report_path)
            print(f"  mapping report: {report_path}")
        print("  dry run: no IMAS entry written" if args.dry_run else f"  wrote: {result.uri}")

    report_path = Path(args.report).expanduser()
    write_batch_report(results, report_path)
    print(f"Batch report: {report_path}")
    uri_map_path = Path(args.uri_map).expanduser()
    write_uri_map(results, uri_map_path)
    print(f"URI map: {uri_map_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
