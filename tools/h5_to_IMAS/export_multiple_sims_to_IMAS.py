#!/usr/bin/env python3
"""Export combined DREAM HDF5 simulations to IMAS with the new mapper path.

Run ``python export_multiple_sims_to_IMAS.py --help`` to see the full workflow,
examples, and command-line options.
"""

from __future__ import annotations

import argparse
import fnmatch
import getpass
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np

import imas_mapper
from dream_h5_reader import read_h5_raw_combined
from imas_mapper import DEFAULT_IDS, build_ids, build_imas_data


DEFAULT_URI_TEMPLATE = "imas:hdf5?user={user};pulse={pulse};run={run};database={database};version={version}"
STAGE_ORDER = {
    "init": 0,
    "injection_start": 10,
    "injection_continuation": 20,
    "injection": 25,
    "cq": 30,
    "unknown": 100,
}


@dataclass
class SimulationResult:
    name: str
    directory: Path
    files: list[Path]
    uri: str
    pulse: int
    run: int
    ids_written: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    time: np.ndarray | None = None
    timings: dict[str, float] = field(default_factory=dict)
    dry_run: bool = False


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


def keep_stage_file(path: Path, include_init: bool = False) -> bool:
    return include_init or classify_stage(path) != "init"


def stage_files_in(directory: Path, stage_glob: str, include_init: bool = False) -> list[Path]:
    return [
        path
        for path in sort_stage_files(directory.glob(stage_glob))
        if path.is_file() and keep_stage_file(path, include_init)
    ]


def discover_simulations(
    root: Path,
    stage_glob: str = "output*.h5",
    include_init: bool = False,
) -> list[tuple[str, Path, list[Path]]]:
    """Find DREAM simulations below root.

    A simulation can be:
    - a direct folder containing matching HDF5 files
    - a case folder with an ``output`` child containing matching HDF5 files
    - one child folder per simulation below ``root``
    """
    simulations: list[tuple[str, Path, list[Path]]] = []

    output_dirs = {
        path.parent
        for path in root.rglob(stage_glob)
        if path.is_file() and path.parent.name == "output"
    }
    for output_dir in sorted(output_dirs, key=lambda path: natural_relative_key(path.parent, root)):
        directory = output_dir.parent
        files = stage_files_in(output_dir, stage_glob, include_init)
        if files:
            simulations.append((simulation_name(root, directory), directory, files))
    if simulations:
        return simulations

    direct_files = stage_files_in(root, stage_glob, include_init)
    if direct_files:
        directory = root.parent if root.name == "output" else root
        return [(simulation_name(root, directory), directory, direct_files)]

    for child in sorted((path for path in root.iterdir() if path.is_dir()), key=lambda path: natural_relative_key(path, root)):
        files = stage_files_in(child, stage_glob, include_init)
        if files:
            simulations.append((simulation_name(root, child), child, files))
    return simulations


def read_folder_filter_file(path: str | None) -> set[str] | None:
    if not path:
        return None
    filter_path = Path(path).expanduser()
    entries: set[str] = set()
    for line in filter_path.read_text(encoding="utf-8").splitlines():
        text = line.split("#", 1)[0].strip()
        if text:
            entries.update(expand_filter_entry(text))
    return entries


def normalize_filter_text(text: str) -> str:
    return Path(text.strip()).as_posix().lstrip("./").rstrip("/")


def expand_filter_entry(text: str) -> set[str]:
    entry = normalize_filter_text(text.split()[0])
    if not entry:
        return set()

    variants = {entry}
    parts = [part for part in entry.split("/") if part]
    if parts and parts[-1] == "output":
        parts = parts[:-1]
        if parts:
            variants.add("/".join(parts))

    for index in range(len(parts)):
        variants.add("/".join(parts[index:]))
    if parts:
        variants.add(parts[-1])
        variants.add("/".join(parts + ["output"]))
    return {variant for variant in variants if variant}


def simulation_filter_keys(root: Path, name: str, directory: Path) -> set[str]:
    keys = {
        normalize_filter_text(name),
        normalize_filter_text(directory.name),
        normalize_filter_text(directory.as_posix()),
        normalize_filter_text((directory / "output").as_posix()),
    }
    try:
        keys.add(normalize_filter_text(directory.relative_to(root).as_posix()))
    except ValueError:
        pass
    try:
        keys.add(normalize_filter_text((directory / "output").relative_to(root).as_posix()))
    except ValueError:
        pass
    try:
        keys.add(normalize_filter_text(directory.relative_to(Path.cwd()).as_posix()))
    except ValueError:
        pass
    try:
        keys.add(normalize_filter_text((directory / "output").relative_to(Path.cwd()).as_posix()))
    except ValueError:
        pass
    return {key for key in keys if key}


def matches_filter_entry(keys: set[str], entries: set[str]) -> bool:
    for entry in entries:
        for key in keys:
            if (
                key == entry
                or key.endswith(f"/{entry}")
                or entry.endswith(f"/{key}")
                or key.startswith(f"{entry}/")
                or f"/{entry}/" in key
                or fnmatch.fnmatch(key, entry)
                or fnmatch.fnmatch(entry, key)
            ):
                return True
    return False


def filter_simulations(
    simulations: list[tuple[str, Path, list[Path]]],
    root: Path,
    include_entries: set[str] | None = None,
    exclude_entries: set[str] | None = None,
) -> list[tuple[str, Path, list[Path]]]:
    filtered = []
    for name, directory, files in simulations:
        keys = simulation_filter_keys(root, name, directory)
        if include_entries is not None and not matches_filter_entry(keys, include_entries):
            continue
        if exclude_entries is not None and matches_filter_entry(keys, exclude_entries):
            continue
        filtered.append((name, directory, files))
    return filtered


def simulation_name(root: Path, directory: Path) -> str:
    try:
        return directory.relative_to(root).as_posix()
    except ValueError:
        return directory.name


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


def export_single_sim(
    name: str,
    directory: Path,
    files: list[Path],
    uri: str,
    pulse: int,
    run: int,
    ids: list[str],
    dd_version: str | None = None,
    dry_run: bool = False,
) -> SimulationResult:
    """Read, combine, map, and optionally write one DREAM simulation."""
    result = SimulationResult(
        name=name,
        directory=directory,
        files=files,
        uri=uri,
        pulse=pulse,
        run=run,
        dry_run=dry_run,
    )

    start = time.perf_counter()
    raw, missing = read_h5_raw_combined(files)
    result.timings["read_combine"] = time.perf_counter() - start
    result.missing = missing
    result.time = np.asarray(raw.get("time"), dtype=float).reshape(-1)

    start = time.perf_counter()
    if dry_run:
        build_imas_data(raw)
        result.ids_written = list(ids)
        result.timings["transform"] = time.perf_counter() - start
        return result

    ids_list = build_ids(raw, selected=ids, dd_version=dd_version)
    result.timings["map_ids"] = time.perf_counter() - start

    start = time.perf_counter()
    imas_mapper.ensure_imas()
    db = imas_mapper.imas.DBEntry(uri, "w")
    try:
        for ids_obj in ids_list:
            db.put(ids_obj)
        result.ids_written = list(ids)
    finally:
        try:
            db.close()
        except Exception:
            pass
    result.timings["write"] = time.perf_counter() - start

    return result


def time_range_text(time_values: np.ndarray | None) -> str:
    if time_values is None or time_values.size == 0:
        return "(empty)"
    return f"{float(time_values[0]):.9g} to {float(time_values[-1]):.9g} s ({time_values.size})"


def write_batch_report(results: list[SimulationResult], path: Path) -> None:
    lines = ["# DREAM IMAS Export Report", ""]
    for result in results:
        lines.extend([
            f"## {result.name}",
            "",
            f"- Directory: `{result.directory}`",
            f"- URI: `{result.uri}`",
            f"- Pulse/run: `{result.pulse}/{result.run}`",
            f"- Stages: {len(result.files)}",
            f"- Time: {time_range_text(result.time)}",
            f"- IDSs: {', '.join(f'`{name}`' for name in result.ids_written) if result.ids_written else '(none)'}",
            f"- Missing aliases: {len(result.missing)}",
            "",
            "| Stage | Type | File |",
            "| --- | --- | --- |",
        ])
        for index, file_path in enumerate(result.files, start=1):
            lines.append(f"| {index} | `{classify_stage(file_path)}` | `{file_path.name}` |")
        lines.append("")
        if result.timings:
            lines.extend(["| Timing Section | Seconds |", "| --- | --- |"])
            for key, seconds in sorted(result.timings.items(), key=lambda item: item[1], reverse=True):
                lines.append(f"| `{key}` | {seconds:.3f} |")
            lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_uri_map(results: list[SimulationResult], path: Path) -> None:
    lines = ["URI\tProvenance folder"]
    lines.extend(f"{result.uri}\t{result.name}" for result in results)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export one or more DREAM HDF5 simulations to IMAS using the new "
            "reader -> mapper -> put() workflow."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""\
What this script does
---------------------
1. Discovers DREAM .h5 stage files below ROOT.
2. Sorts each simulation's stages in a deterministic DREAM order:
   init, injection_start, injection_continuation, injection, cq, then unknown.
3. Reads the selected HDF5 files with dream_h5_reader.py.
4. Combines stage dictionaries into one raw simulation dictionary.
5. Builds IMAS-ready static/dynamic data and IDS objects with imas_mapper.py.
6. Writes each selected IDS with db.put(ids), not put_slice().
7. Writes a Markdown report and a TSV URI/provenance map.

Simulation discovery
--------------------
ROOT may be:
  - a directory directly containing files matching --stage-glob
  - a case directory containing an output/ child with matching files
  - a parent directory containing one subdirectory per simulation

By default, files classified as init/initial are excluded from the time-appended
stage list. Use --include-init only if you deliberately want init files to be
combined as time stages.

Folder include/exclude files
----------------------------
--include-folders-file and --exclude-folders-file read plain text files with one
simulation folder per line. Blank lines and text after # are ignored.

Entries may be:
  - the discovered simulation name relative to ROOT
  - the simulation folder basename
  - a path relative to ROOT
  - an absolute folder path
  - the same folder path with a trailing output component

If --include-folders-file is given, only matching simulations are exported.
If --exclude-folders-file is given, matching simulations are removed after the
include filter is applied. Matching is suffix-aware, so the same filter file can
be used when ROOT is a high-level parent, Baseline, St4, or St4/output.

Default IDS list
----------------
{', '.join(DEFAULT_IDS)}

URI template fields
-------------------
The --uri-template string can use:
  {{sim}}       simulation name relative to ROOT
  {{case}}      same as {{sim}}
  {{sim_dir}}   absolute simulation directory
  {{index}}     1-based simulation index in this batch
  {{pulse}}     assigned pulse number
  {{run}}       assigned run number
  {{user}}      --imas-user
  {{database}}  --database
  {{version}}   --version

Examples
--------
Dry-run one simulation folder and only validate reading/combining/transforms:
  python export_multiple_sims_to_IMAS.py /path/to/case/output --dry-run

Export all simulations below a parent directory with default IDSs:
  python export_multiple_sims_to_IMAS.py /path/to/simulations

Export only selected IDSs:
  python export_multiple_sims_to_IMAS.py /path/to/simulations --ids summary plasma_profiles spi

Use a custom stage glob and pulse range:
  python export_multiple_sims_to_IMAS.py /path/to/simulations --stage-glob '*.h5' --pulse-start 130000 --pulse-end 130100

Use a custom IMAS URI template:
  python export_multiple_sims_to_IMAS.py /path/to/simulations \\
    --uri-template 'imas:hdf5?user={{user}};pulse={{pulse}};run={{run}};database={{database}};version={{version}}'

Outputs
-------
--report writes a human-readable Markdown report with stage order, time range,
IDS list, missing aliases, and timings.

--uri-map writes a TSV file with the IMAS URI and provenance folder for each
exported simulation.
""",
    )
    parser.add_argument(
        "root",
        help=(
            "Simulation folder, output folder, or parent folder containing multiple "
            "simulation folders."
        ),
    )
    parser.add_argument(
        "--stage-glob",
        default="output*.h5",
        help=(
            "Glob used to find DREAM output stages inside candidate directories. "
            "Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--include-init",
        action="store_true",
        help=(
            "Include init/initial files in the time-appended stage list. By default "
            "init files are skipped so they do not create an artificial time stage."
        ),
    )
    parser.add_argument(
        "--include-folders-file",
        default=None,
        help=(
            "Text file listing simulation folders to include. One entry per line; "
            "blank lines and # comments are ignored. Entries may be discovered "
            "simulation names, basenames, paths relative to ROOT, absolute paths, "
            "or paths ending in output."
        ),
    )
    parser.add_argument(
        "--exclude-folders-file",
        default=None,
        help=(
            "Text file listing simulation folders to exclude after the include "
            "filter is applied. Format is the same as --include-folders-file."
        ),
    )
    parser.add_argument(
        "--uri-template",
        default=DEFAULT_URI_TEMPLATE,
        help=(
            "IMAS URI template. Available fields: {sim}, {case}, {sim_dir}, {index}, "
            "{pulse}, {run}, {user}, {database}, {version}."
        ),
    )
    parser.add_argument(
        "--imas-user",
        default=getpass.getuser(),
        help="IMAS user field used by the default URI template. Default: %(default)s.",
    )
    parser.add_argument(
        "--database",
        default="ITER_DISRUPTIONS_DREAM",
        help="IMAS database field used by the default URI template. Default: %(default)s.",
    )
    parser.add_argument(
        "--version",
        default="4",
        help="IMAS version field used by the default URI template. Default: %(default)s.",
    )
    parser.add_argument(
        "--pulse-start",
        type=int,
        default=120000,
        help=(
            "First pulse number assigned to discovered simulations. Each simulation "
            "gets the next pulse number. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--pulse-end",
        type=int,
        default=129999,
        help=(
            "Maximum pulse number allowed for this export. The script stops before "
            "writing if the discovered simulations exceed this range. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--run",
        type=int,
        default=1,
        help="Run number assigned to each exported simulation. Default: %(default)s.",
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        default=DEFAULT_IDS,
        choices=DEFAULT_IDS,
        help=(
            "IDSs to create/write. Choices: %(choices)s. Default: all currently "
            "implemented IDSs."
        ),
    )
    parser.add_argument(
        "--dd-version",
        default=None,
        help=(
            "Optional IMAS Data Dictionary version passed to imas.IDSFactory. "
            "Default: use the environment/default IMAS version."
        ),
    )
    parser.add_argument(
        "--report",
        default="dream_multi_stage_export_report2.md",
        help=(
            "Markdown batch report path. Contains discovered stages, time ranges, "
            "missing aliases, selected IDSs, URI, and timings. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--uri-map",
        default="dream_simulation_uris2.tsv",
        help=(
            "TSV file listing each written URI and the simulation/provenance folder. "
            "Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Read, combine, and transform data, but do not instantiate IDS objects "
            "and do not write to IMAS. Useful for testing discovery and raw data coverage."
        ),
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Print coarse timing sections for each simulation while running.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        print(f"ERROR: root path does not exist: {root}", file=sys.stderr)
        return 2

    simulations = discover_simulations(root, args.stage_glob, args.include_init)
    if not simulations:
        print(f"ERROR: no stage files matching {args.stage_glob!r} found under {root}", file=sys.stderr)
        return 2
    include_entries = read_folder_filter_file(args.include_folders_file)
    exclude_entries = read_folder_filter_file(args.exclude_folders_file)
    if include_entries is not None or exclude_entries is not None:
        discovered_count = len(simulations)
        discovered_sample = simulations[:5]
        simulations = filter_simulations(simulations, root, include_entries, exclude_entries)
        if not simulations:
            sample_names = [
                "/".join(sorted(simulation_filter_keys(root, name, directory))[:3])
                for name, directory, _ in discovered_sample
            ]
            print(
                "ERROR: folder filters removed all discovered simulations "
                f"({discovered_count} before filtering).",
                file=sys.stderr,
            )
            print(
                "       Check that entries in --include-folders-file/--exclude-folders-file "
                "match simulation folders, not individual HDF5 files.",
                file=sys.stderr,
            )
            if include_entries:
                print(f"       Include entries sample: {', '.join(sorted(include_entries)[:5])}", file=sys.stderr)
            if exclude_entries:
                print(f"       Exclude entries sample: {', '.join(sorted(exclude_entries)[:5])}", file=sys.stderr)
            if sample_names:
                print(f"       Discovered simulation keys sample: {', '.join(sample_names)}", file=sys.stderr)
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
        uri, pulse, run = simulation_uri(name, directory, args, index)
        print(f"[{index}/{len(simulations)}] {name}: {len(files)} stage file(s)")
        print("  stage order: " + " -> ".join(f"{classify_stage(path)}:{path.name}" for path in files))
        result = export_single_sim(
            name=name,
            directory=directory,
            files=files,
            uri=uri,
            pulse=pulse,
            run=run,
            ids=list(args.ids),
            dd_version=args.dd_version,
            dry_run=args.dry_run,
        )
        results.append(result)
        if args.profile:
            for key, seconds in sorted(result.timings.items(), key=lambda item: item[1], reverse=True):
                print(f"  {key}: {seconds:.3f}s")
        print("  dry run: no IMAS entry written" if args.dry_run else f"  wrote: {result.uri}")

    write_batch_report(results, Path(args.report).expanduser())
    write_uri_map(results, Path(args.uri_map).expanduser())
    print(f"Batch report: {args.report}")
    print(f"URI map: {args.uri_map}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
