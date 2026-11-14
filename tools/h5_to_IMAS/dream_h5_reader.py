"""Raw DREAM HDF5 reader helpers.

This module owns the DREAM-native dictionary layer:

- read one DREAM HDF5 file into an aliases-only raw dictionary
- combine several raw dictionaries for one simulation

"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

try:
    import h5py
except Exception as exc:  # pragma: no cover - depends on user environment
    h5py = None  # type: ignore[assignment]
    H5PY_IMPORT_ERROR = exc
else:
    H5PY_IMPORT_ERROR = None


H5_ALIAS_PATHS: dict[str, str] = {
    # Common grids and geometry.
    "time": "/grid/t",
    "r": "/grid/r",
    "dr": "/grid/dr",
    "R0": "/grid/R0",
    "R0_eq": "/grid/eq/R0",
    "a": "/grid/a",
    "VpVol": "/grid/VpVol",
    "GR0": "/grid/geometry/GR0",
    "FSA_R02OverR2": "/grid/geometry/FSA_R02OverR2",
    "Bmin": "/grid/geometry/Bmin",
    "Bmax": "/grid/geometry/Bmax",
    "FSA_BOverBmin": "/grid/geometry/FSA_BOverBmin",
    "FSA_BOverBmin2": "/grid/geometry/FSA_BOverBmin2",
    "effectivePassingFraction": "/grid/geometry/effectivePassingFraction",
    "xi0TrappedBoundary": "/grid/geometry/xi0TrappedBoundary",
    "toroidalFlux": "/grid/geometry/toroidalFlux",
    # Equilibrium geometry.
    "Z0_eq": "/grid/eq/Z0",
    "theta_eq": "/grid/eq/theta",
    "RMinusR0": "/grid/eq/RMinusR0",
    "ZMinusZ0": "/grid/eq/ZMinusZ0",
    "RMinusR0_f": "/grid/eq/RMinusR0_f",
    "ZMinusZ0_f": "/grid/eq/ZMinusZ0_f",
    # Equation-system quantities.
    "T_cold": "/eqsys/T_cold",
    "n_tot": "/eqsys/n_tot",
    "n_cold": "/eqsys/n_cold",
    "n_i": "/eqsys/n_i",
    "n_re": "/eqsys/n_re",
    "j_tot": "/eqsys/j_tot",
    "j_ohm": "/eqsys/j_ohm",
    "j_re": "/eqsys/j_re",
    "E_field": "/eqsys/E_field",
    "I_p": "/eqsys/I_p",
    "psi_p": "/eqsys/psi_p",
    "W_cold": "/eqsys/W_cold",
    "W_i": "/eqsys/W_i",
    "x_p": "/eqsys/x_p",
    "v_p": "/eqsys/v_p",
    "Y_p": "/eqsys/Y_p",
    # Other DREAM output groups.
    "conductivity": "/other/fluid/conductivity",
    "Zeff": "/other/fluid/Zeff",
    "runawayRate": "/other/fluid/runawayRate",
    "gammaDreicer": "/other/fluid/gammaDreicer",
    "gammaHottail": "/other/fluid/gammaHottail",
    "gammaTritium": "/other/fluid/gammaTritium",
    "gammaCompton": "/other/fluid/gammaCompton",
    "pCrit": "/other/fluid/pCrit",
    "pCritHottail": "/other/fluid/pCritHottail",
    "EDreic": "/other/fluid/EDreic",
    "Ectot": "/other/fluid/Ectot",
    "Tcold_radiation": "/other/fluid/Tcold_radiation",
    # Ion and SPI settings.
    "ion_Z": "/settings/eqsys/n_i/Z",
    "ion_isotopes": "/settings/eqsys/n_i/isotopes",
    "ion_names": "/settings/eqsys/n_i/names",
    "ion_hydrogennames": "/settings/eqsys/n_i/hydrogennames",
    "ion_tritiumnames": "/settings/eqsys/n_i/tritiumnames",
    "SPIMolarFraction": "/settings/eqsys/n_i/SPIMolarFraction",
    "spi_Ninj": "/settings/eqsys/spi/init/Ninj",
    "spi_ZsDrift": "/settings/eqsys/spi/ZsDrift",
    "spi_isotopesDrift": "/settings/eqsys/spi/isotopesDrift",
    # Code provenance.
    "code_commit": "/code/commit",
    "code_changes": "/code/changes",
    "code_refspec": "/code/refspec",
    # Older DREAM metadata fallback paths used by the current mapper.
    "ionmeta_Z": "/ionmeta/Z",
    "ionmeta_names": "/ionmeta/names",
}

TIME_DEPENDENT_ALIASES: set[str] = {
    "T_cold",
    "n_tot",
    "n_cold",
    "n_i",
    "n_re",
    "j_tot",
    "j_ohm",
    "j_re",
    "E_field",
    "I_p",
    "psi_p",
    "W_cold",
    "W_i",
    "x_p",
    "v_p",
    "Y_p",
    "conductivity",
    "Zeff",
    "runawayRate",
    "gammaDreicer",
    "gammaHottail",
    "gammaTritium",
    "gammaCompton",
    "pCrit",
    "pCritHottail",
    "EDreic",
    "Ectot",
    "Tcold_radiation",
}

TIME_TOLERANCE = 1e-12

__all__ = [
    "H5_ALIAS_PATHS",
    "TIME_DEPENDENT_ALIASES",
    "combine_raw_dicts",
    "decode_h5_value",
    "read_h5_raw",
    "read_h5_raw_combined",
]


def decode_h5_value(value: Any) -> Any:
    """Decode HDF5 string-like values while preserving numeric arrays."""
    arr = np.asarray(value)

    if arr.dtype.kind == "S":
        if arr.ndim == 0:
            return arr.item().decode(errors="ignore").rstrip("\x00")
        chars = [item.decode(errors="ignore") for item in arr.reshape(-1)]
        return "".join(chars).rstrip("\x00")

    if arr.dtype.kind == "U":
        if arr.ndim == 0:
            return str(arr.item()).rstrip("\x00")
        return "".join(str(item) for item in arr.reshape(-1)).rstrip("\x00")

    if arr.dtype.kind == "O" and arr.size > 0:
        flat = arr.reshape(-1)
        if all(isinstance(item, (bytes, str)) for item in flat):
            text = "".join(
                item.decode(errors="ignore") if isinstance(item, bytes) else str(item)
                for item in flat
            )
            return text.rstrip("\x00")

    return arr


def read_h5_raw(path: str | Path) -> tuple[dict[str, Any], list[str]]:
    """Read one DREAM HDF5 file into an aliases-only raw dictionary.

    Missing datasets are represented as ``None`` in the returned dictionary and
    as ``"<alias>: <hdf5 path>"`` entries in the missing list.
    """
    h5_path = Path(path).expanduser()
    if not h5_path.exists():
        raise FileNotFoundError(f"DREAM HDF5 file does not exist: {h5_path}")

    if h5py is None:
        raise RuntimeError(
            "Could not import h5py. Install h5py or load the DREAM Python environment. "
            f"Original import error: {H5PY_IMPORT_ERROR!r}"
        )

    data: dict[str, Any] = {}
    missing: list[str] = []

    with h5py.File(h5_path, "r") as h5:
        for alias, dataset_path in H5_ALIAS_PATHS.items():
            if dataset_path not in h5:
                data[alias] = None
                missing.append(f"{alias}: {dataset_path}")
                continue
            data[alias] = decode_h5_value(h5[dataset_path][()])

    return data, missing


def read_h5_raw_combined(paths: list[str | Path]) -> tuple[dict[str, Any], list[str]]:
    """Read several DREAM HDF5 files and combine them into one raw dictionary."""
    stage_data = []
    missing = []

    for path in paths:
        data, stage_missing = read_h5_raw(path)
        stage_data.append(data)
        missing.extend(f"{Path(path).name}: {item}" for item in stage_missing)

    return combine_raw_dicts(stage_data), missing


def combine_raw_dicts(stage_data: list[dict[str, Any]]) -> dict[str, Any]:
    """Combine raw dictionaries from several DREAM stages.

    Time-dependent aliases are concatenated along axis 0. Static aliases are
    kept from the first stage and must match in all later stages where present.
    The returned ``time`` array is shifted stage-by-stage so it is strictly
    increasing across the combined simulation.
    """
    if not stage_data:
        raise ValueError("Cannot combine DREAM stages: no input dictionaries were provided.")

    stage_times = [_time_array(stage, index) for index, stage in enumerate(stage_data)]
    corrected_times = _correct_stage_times(stage_times)

    combined: dict[str, Any] = {}
    for alias in H5_ALIAS_PATHS:
        if alias == "time":
            combined[alias] = np.concatenate(corrected_times) if corrected_times else np.array([])
        elif alias in TIME_DEPENDENT_ALIASES:
            combined[alias] = _combine_time_dependent(alias, stage_data, stage_times)
        else:
            combined[alias] = _combine_static(alias, stage_data)

    return combined

# Functions to correctly combine time arrays
def _time_array(stage: dict[str, Any], index: int) -> np.ndarray:
    time = stage.get("time")
    if time is None:
        raise ValueError(f"Cannot combine DREAM stages: stage {index} is missing the time array.")
    arr = np.asarray(time, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"Cannot combine DREAM stages: stage {index} has an empty time array.")
    return arr


def _correct_stage_times(stage_times: list[np.ndarray]) -> list[np.ndarray]:
    corrected = []
    previous_end: float | None = None

    for time in stage_times:
        shifted = np.array(time, dtype=float, copy=True)
        if previous_end is not None and shifted[0] <= previous_end + TIME_TOLERANCE:
            shifted += previous_end - shifted[0] + _positive_time_step(shifted)
        corrected.append(shifted)
        previous_end = float(shifted[-1])

    return corrected


def _positive_time_step(time: np.ndarray) -> float:
    diffs = np.diff(time)
    positive = diffs[diffs > TIME_TOLERANCE]
    if positive.size > 0:
        return float(positive[0])
    return TIME_TOLERANCE

# Functions to combine time-dependent and static aliases
def _combine_time_dependent(
    alias: str,
    stage_data: list[dict[str, Any]],
    stage_times: list[np.ndarray],
) -> Any:
    values = [stage.get(alias) for stage in stage_data]
    if all(value is None for value in values):
        return None
    if any(value is None for value in values):
        raise ValueError(f"{alias}: time-dependent data is missing in at least one stage.")

    arrays = []
    for index, (value, time) in enumerate(zip(values, stage_times)):
        arr = np.asarray(value)
        if arr.ndim == 0:
            raise ValueError(f"{alias}: stage {index} is scalar, expected a time-dependent array.")
        valid_time_lengths = {time.size}
        if time.size > 1:
            valid_time_lengths.add(time.size - 1)
        if arr.shape[0] not in valid_time_lengths:
            raise ValueError(
                f"{alias}: stage {index} first dimension is {arr.shape[0]}, "
                f"expected {sorted(valid_time_lengths)} from its time array."
            )
        arrays.append(arr)

    try:
        return np.concatenate(arrays, axis=0)
    except ValueError as exc:
        raise ValueError(f"{alias}: could not concatenate stage arrays: {exc}") from exc


def _combine_static(alias: str, stage_data: list[dict[str, Any]]) -> Any:
    values = [stage.get(alias) for stage in stage_data if stage.get(alias) is not None]
    if not values:
        return None

    reference = values[0]
    for value in values[1:]:
        if not _values_equal(reference, value):
            raise ValueError(f"{alias}: static data differs between DREAM stages.")

    return reference


def _values_equal(left: Any, right: Any) -> bool:
    left_arr = np.asarray(left)
    right_arr = np.asarray(right)

    if left_arr.shape != right_arr.shape:
        return False
    if left_arr.dtype.kind in {"S", "U", "O"} or right_arr.dtype.kind in {"S", "U", "O"}:
        return bool(np.array_equal(left_arr, right_arr))
    return bool(np.array_equal(left_arr, right_arr, equal_nan=True))
