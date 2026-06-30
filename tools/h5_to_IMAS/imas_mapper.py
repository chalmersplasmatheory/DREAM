"""Transform combined DREAM dictionaries and map them to IMAS IDS objects.

This module starts after dream_h5_reader has already read and combined the
DREAM HDF5 files. It does not read HDF5 files directly.

Filled IDS fields and DREAM origins
-----------------------------------

Shared grid and metadata helpers:

*.vacuum_toroidal_field.r0 | /grid/R0.
*.vacuum_toroidal_field.b0 | time trace filled with B0 derived from PHI_SIGN * /grid/geometry/GR0[-1].
*.profiles_1d.grid.rho_tor | derived from sqrt(/grid/geometry/toroidalFlux / (pi * abs(B0))).
*.profiles_1d.grid.rho_tor_norm | normalized rho_tor.
*.profiles_1d.grid.psi | PSI_COCOS * /eqsys/psi_p.
*.profiles_1d.grid.psi_magnetic_axis | first radial value of mapped psi.
*.profiles_1d.grid.psi_boundary | last radial value of mapped psi.
*.profiles_1d.grid.rho_pol_norm | normalized mapped psi.
*.code.name | constant DREAM.
*.code.description | constant DREAM description.
*.code.commit | /code/commit.
*.code.repository | constant DREAM public repository URL.
*.code.version | /code/refspec.
*.code.parameters | XML fragment built from /code/commit, /code/changes, /code/datetime_commit, /code/datetime_simulation, and /code/refspec.

summary:

summary.time | corrected combined /grid/t.
summary.global_quantities.r0.value | /grid/R0.
summary.global_quantities.b0.value | time trace filled with derived B0.
summary.global_quantities.ip.value | PHI_SIGN * /eqsys/I_p.
summary.global_quantities.current_ohm.value | radial integral of /eqsys/j_ohm with DREAM area weights.
summary.global_quantities.energy_electrons_thermal.value | volume integral of /eqsys/W_cold.
summary.global_quantities.energy_ion_total_thermal.value | volume integral of /eqsys/W_i.
summary.global_quantities.energy_thermal.value | sum of electron and ion thermal energies when available.
summary.volume_average.t_e.value | volume average of /eqsys/T_cold.
summary.volume_average.n_e.value | volume average of /eqsys/n_tot; falls back to /eqsys/n_cold.
summary.volume_average.zeff.value | volume average of /other/fluid/Zeff.
summary.local.magnetic_axis.t_e.value | first radial sample of /eqsys/T_cold.
summary.local.magnetic_axis.n_e.value | first radial sample of electron density.
summary.local.magnetic_axis.zeff.value | first radial sample of /other/fluid/Zeff.
summary.local.magnetic_axis.e_field_parallel.value | first radial sample of /eqsys/E_field.
summary.local.separatrix.t_e.value | last radial sample of /eqsys/T_cold.
summary.local.separatrix.n_e.value | last radial sample of electron density.
summary.local.separatrix.zeff.value | last radial sample of /other/fluid/Zeff.
summary.local.separatrix.e_field_parallel.value | last radial sample of /eqsys/E_field.
summary.runaways.current.value | radial integral of PHI_SIGN * /eqsys/j_re with DREAM area weights.
summary.runaways.current_phi_max.value | maximum absolute runaway current.
summary.runaways.particles.value | volume integral of /eqsys/n_re.

runaway_electrons:

runaway_electrons.time | corrected combined /grid/t.
runaway_electrons.profiles_1d[*].density | /eqsys/n_re.
runaway_electrons.profiles_1d[*].current_density |  /eqsys/j_re.
runaway_electrons.profiles_1d[*].ddensity_dt_total | /other/fluid/runawayRate.
runaway_electrons.profiles_1d[*].ddensity_dt_dreicer | /other/fluid/gammaDreicer.
runaway_electrons.profiles_1d[*].ddensity_dt_hot_tail | /other/fluid/gammaHottail.
runaway_electrons.profiles_1d[*].ddensity_dt_tritium | /other/fluid/gammaTritium.
runaway_electrons.profiles_1d[*].ddensity_dt_compton | /other/fluid/gammaCompton.
runaway_electrons.profiles_1d[*].momentum_critical_avalanche | M_ELECTRON * C_LIGHT * /other/fluid/pCrit.
runaway_electrons.profiles_1d[*].momentum_critical_hot_tail | M_ELECTRON * C_LIGHT * /other/fluid/pCritHottail.
runaway_electrons.profiles_1d[*].e_field_dreicer | /other/fluid/EDreic.
runaway_electrons.profiles_1d[*].e_field_critical | /other/fluid/Ectot.
runaway_electrons.global_quantities.current_phi | radial integral of PHI_SIGN * /eqsys/j_re with DREAM area weights.

plasma_profiles:

plasma_profiles.time | corrected combined /grid/t.
plasma_profiles.profiles_1d[*].electrons.temperature | /eqsys/T_cold.
plasma_profiles.profiles_1d[*].electrons.density | /eqsys/n_tot; falls back to /eqsys/n_cold.
plasma_profiles.profiles_1d[*].electrons.density_thermal | /eqsys/n_cold.
plasma_profiles.profiles_1d[*].conductivity_parallel | /other/fluid/conductivity.
plasma_profiles.profiles_1d[*].e_field.parallel | /eqsys/E_field.
plasma_profiles.profiles_1d[*].zeff | /other/fluid/Zeff.
plasma_profiles.profiles_1d[*].j_total | /eqsys/j_tot.
plasma_profiles.profiles_1d[*].j_ohmic | /eqsys/j_ohm.
plasma_profiles.profiles_1d[*].ion[*].neutral_index | DREAM ion species index.
plasma_profiles.profiles_1d[*].neutral[*].ion_index | DREAM ion species index.
plasma_profiles.profiles_1d[*].neutral[*].element[0].z_n | /settings/eqsys/n_i/Z or /ionmeta/Z.
plasma_profiles.profiles_1d[*].neutral[*].element[0].a | /settings/eqsys/n_i/isotopes with H/D/T fallback rules.
plasma_profiles.profiles_1d[*].neutral[*].density | /eqsys/n_i neutral charge state Z0=0.
plasma_profiles.profiles_1d[*].neutral[*].state[0].density | /eqsys/n_i neutral charge state Z0=0.
plasma_profiles.profiles_1d[*].ion[*].name | /settings/eqsys/n_i/names or /ionmeta/names reconciled with Z.
plasma_profiles.profiles_1d[*].ion[*].element[0].z_n | /settings/eqsys/n_i/Z or /ionmeta/Z.
plasma_profiles.profiles_1d[*].ion[*].element[0].a | /settings/eqsys/n_i/isotopes with H/D/T fallback rules.
plasma_profiles.profiles_1d[*].ion[*].density | sum of /eqsys/n_i charged states Z0=1..Z.
plasma_profiles.profiles_1d[*].ion[*].state[*].name | charge-state label "<species> Z0=<charge>".
plasma_profiles.profiles_1d[*].ion[*].state[*].density | /eqsys/n_i charged state Z0=1..Z.

equilibrium:

equilibrium.time | corrected combined /grid/t.
equilibrium.time_slice[*].global_quantities.ip | PHI_SIGN * /eqsys/I_p.
equilibrium.time_slice[*].profiles_1d.rho_tor | derived from toroidal flux and B0.
equilibrium.time_slice[*].profiles_1d.rho_tor_norm | normalized rho_tor.
equilibrium.time_slice[*].profiles_1d.phi | /grid/geometry/toroidalFlux.
equilibrium.time_slice[*].profiles_1d.r_outboard | /grid/R0 + /grid/r.
equilibrium.time_slice[*].profiles_1d.b_field_min | PHI_SIGN * /grid/geometry/Bmin.
equilibrium.time_slice[*].profiles_1d.b_field_max | PHI_SIGN * /grid/geometry/Bmax.
equilibrium.time_slice[*].profiles_1d.b_field_average | b_field_min * /grid/geometry/FSA_BOverBmin.
equilibrium.time_slice[*].profiles_1d.trapped_fraction | 1 - /grid/geometry/effectivePassingFraction, clipped to [0, 1].
equilibrium.time_slice[*].profiles_1d.gm1 | /grid/geometry/FSA_R02OverR2 / R0**2.
equilibrium.time_slice[*].profiles_1d.gm5 | b_field_min**2 * /grid/geometry/FSA_BOverBmin2.
equilibrium.time_slice[*].profiles_1d.j_parallel | /eqsys/j_tot.
equilibrium.time_slice[*].profiles_1d.psi | PSI_COCOS * /eqsys/psi_p.
equilibrium.time_slice[*].profiles_1d.psi_norm | normalized mapped psi.
equilibrium.time_slice[*].global_quantities.psi_magnetic_axis | first radial value of mapped psi.
equilibrium.time_slice[*].global_quantities.psi_boundary | last radial value of mapped psi.
equilibrium.time_slice[*].profiles_1d.dvolume_dpsi | derived from dVdr = /grid/VpVol * R0, /grid/dr, and mapped psi cell edges.
equilibrium.time_slice[*].boundary.outline.r | R0 + /grid/eq/RMinusR0_f[:,-1]; falls back to RMinusR0.
equilibrium.time_slice[*].boundary.outline.z | Z0 + /grid/eq/ZMinusZ0_f[:,-1]; falls back to ZMinusZ0.
equilibrium.time_slice[*].boundary.geometric_axis.r | midpoint of boundary r extent.
equilibrium.time_slice[*].boundary.geometric_axis.z | midpoint of boundary z extent.
equilibrium.time_slice[*].profiles_2d[0].type.name | constant total.
equilibrium.time_slice[*].profiles_2d[0].type.index | constant 0.
equilibrium.time_slice[*].profiles_2d[0].type.description | constant Total fields.
equilibrium.time_slice[*].profiles_2d[0].grid_type.name | constant inverse_psi_polar.
equilibrium.time_slice[*].profiles_2d[0].grid_type.index | constant 13.
equilibrium.time_slice[*].profiles_2d[0].grid_type.description | constant grid description.
equilibrium.time_slice[*].profiles_2d[0].grid.dim1 | mapped psi; cell-edge psi when DREAM R/Z geometry is on edges.
equilibrium.time_slice[*].profiles_2d[0].grid.dim2 | /grid/eq/theta; falls back to uniform 0..2*pi.
equilibrium.time_slice[*].profiles_2d[0].r | R0 + /grid/eq/RMinusR0_f or RMinusR0.
equilibrium.time_slice[*].profiles_2d[0].z | Z0 + /grid/eq/ZMinusZ0_f or ZMinusZ0.
equilibrium.time_slice[*].profiles_2d[0].psi | mapped psi repeated on the static DREAM R/Z flux-surface geometry.

radiation:

radiation.time | time axis matched to /other/fluid/Tcold_radiation; uses time[1:] for nt-1 DREAM output.
radiation.process[0].identifier.index | constant -1.
radiation.process[0].identifier.name | constant dream_total_electron_radiation.
radiation.process[0].identifier.description | constant DREAM total radiation description.
radiation.process[0].global_quantities[*].inside_vessel.power_electrons | volume integral of /other/fluid/Tcold_radiation.
radiation.process[0].profiles_1d[*].grid.rho_tor_norm | normalized rho_tor.
radiation.process[0].profiles_1d[*].grid.rho_tor | derived rho_tor.
radiation.process[0].profiles_1d[*].electrons.emissivity | /other/fluid/Tcold_radiation.
radiation.process[0].profiles_1d[*].electrons.power_inside | cumulative volume integral of /other/fluid/Tcold_radiation.

spi:

spi.time | corrected combined /grid/t.
spi.injector[*].name | DREAM_SPI or numbered injector name.
spi.injector[*].description | constant DREAM SPI description.
spi.injector[*].pellet.core.atoms_n | group total from /settings/eqsys/spi/init/Ninj and SPIMolarFraction.
spi.injector[*].pellet.core.species[*].name | injected species label from /settings/eqsys/n_i/names or spi/ZsDrift.
spi.injector[*].pellet.core.species[*].z_n | /settings/eqsys/spi/ZsDrift or injected ion Z.
spi.injector[*].pellet.core.species[*].a | /settings/eqsys/spi/isotopesDrift or injected ion isotope.
spi.injector[*].pellet.core.species[*].density | species atoms divided by initial group pellet volume from /eqsys/Y_p or /settings/eqsys/spi/init/rp.
spi.injector[*].velocity_mass_centre_fragments_r | volume-weighted mean cylindrical radial velocity from initial /eqsys/x_p and /eqsys/v_p; falls back to unweighted mean if volumes are missing.
spi.injector[*].velocity_mass_centre_fragments_phi | volume-weighted mean cylindrical toroidal velocity.
spi.injector[*].velocity_mass_centre_fragments_z | volume-weighted mean vertical velocity.
spi.injector[*].fragment[*].position.r | cylindrical R converted from /eqsys/x_p; falls back to repeated spi/init/xp.
spi.injector[*].fragment[*].position.phi | cylindrical phi converted from SPI position.
spi.injector[*].fragment[*].position.z | vertical coordinate from SPI position.
spi.injector[*].fragment[*].velocity_r | cylindrical radial velocity from /eqsys/x_p and /eqsys/v_p; if position is absent, DREAM local radial velocity is used.
spi.injector[*].fragment[*].velocity_phi | cylindrical toroidal velocity; if position is absent, set to zero.
spi.injector[*].fragment[*].velocity_z | vertical velocity from /eqsys/v_p.
spi.injector[*].fragment[*].volume | derived from radius = /eqsys/Y_p**(3/5) or spi/init/rp and volume = 4*pi*radius**3/3.

"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

try:
    import imas
except Exception as exc:  # pragma: no cover - depends on user environment
    imas = None  # type: ignore[assignment]
    IMAS_IMPORT_ERROR = exc
else:
    IMAS_IMPORT_ERROR = None


PSI_COCOS = 2.0 * np.pi
PHI_SIGN = -1.0
C_LIGHT = 299792458.0
M_ELECTRON = 9.10938356e-31
P_NORM = M_ELECTRON * C_LIGHT

DEFAULT_IDS = [
    "summary",
    "runaway_electrons",
    "plasma_profiles",
    "equilibrium",
    "radiation",
    "spi",
]

__all__ = [
    "DEFAULT_IDS",
    "build_ids",
    "build_imas_data",
    "build_static_data",
    "build_dynamic_data",
    "map_equilibrium",
    "map_plasma_profiles",
    "map_radiation",
    "map_runaway_electrons",
    "map_spi",
    "map_summary",
]


@dataclass
class FluxSurface2DContext:
    psi_on_edges: bool
    theta: np.ndarray
    r_2d: np.ndarray
    z_2d: np.ndarray


@dataclass
class SpiComponent:
    name: str
    z: int
    isotope: int | None
    atoms_n: float | None


@dataclass
class SpiInjectorGroup:
    shard_indices: np.ndarray
    fractions: np.ndarray | None


@dataclass
class IonProfileContext:
    n_i: np.ndarray
    z_list: list[int]
    labels: list[str]
    a_list: list[int | None]


def ensure_imas() -> None:
    if imas is None:
        raise RuntimeError(
            "Could not import imas. Install IMAS-Python / load your IMAS environment. "
            f"Original import error: {IMAS_IMPORT_ERROR!r}"
        )


def make_factory(dd_version: str | None = None) -> Any:
    ensure_imas()
    if dd_version:
        return imas.IDSFactory(version=dd_version)
    return imas.IDSFactory()


def make_ids(factory: Any, name: str) -> Any | None:
    if not hasattr(factory, name):
        return None
    ids = getattr(factory, name)()
    set_homogeneous_time(ids)
    try:
        ids.ids_properties.comment = "Generated from combined DREAM HDF5 output."
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


def resize_child_aos(parent: Any, child_name: str, n: int) -> bool:
    if not hasattr(parent, child_name):
        return False
    return resize_aos(getattr(parent, child_name), n)


def fill_ids_field(root: Any, path: str, value: Any) -> bool:
    """Set an IMAS field when both the value and target node exist."""
    if value is None:
        return False

    parts = [part for part in path.split("/") if part]
    if not parts:
        return False

    obj = root
    try:
        for part in parts[:-1]:
            obj = obj[int(part)] if part.isdigit() else getattr(obj, part)
        leaf = parts[-1]
        if not hasattr(obj, leaf):
            return False
        setattr(obj, leaf, sanitize_for_imas(value))
        return True
    except Exception:
        return False


def sanitize_for_imas(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"S", "U", "O"}:
            return value.astype(str)
        return np.asarray(value)
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def build_imas_data(raw: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Convert combined DREAM-native data into static/dynamic IMAS-ready data."""
    static = build_static_data(raw)
    dynamic = build_dynamic_data(raw, static)
    return static, dynamic


def build_static_data(raw: dict[str, Any]) -> dict[str, Any]:
    time = flatten_1d(raw.get("time"))
    r = flatten_1d(raw.get("r"))
    dr = flatten_1d(raw.get("dr"))

    if time is None or time.size == 0:
        raise ValueError("Combined DREAM data does not contain a valid time array.")
    if r is None or r.size == 0:
        raise ValueError("Combined DREAM data does not contain a valid radial grid array.")

    r0 = scalar_first(raw.get("R0"), scalar_first(raw.get("R0_eq")))
    a_minor = scalar_first(raw.get("a"))

    gr0 = scale_optional(flatten_1d(raw.get("GR0")), PHI_SIGN)
    b0 = float(gr0[-1]) if gr0 is not None and gr0.size > 0 else None

    phi_tor = flatten_1d(raw.get("toroidalFlux"))
    rho_tor = None
    rho_tor_norm = None
    if phi_tor is not None and b0 is not None and abs(b0) > 0:
        rho_tor = np.sqrt(phi_tor / (np.pi * abs(b0)))
        rho_tor_norm = normalized_radius(rho_tor)

    vpvol = flatten_1d(raw.get("VpVol"))
    r2inv = flatten_1d(raw.get("FSA_R02OverR2"))
    bmin = scale_optional(flatten_1d(raw.get("Bmin")), PHI_SIGN)
    bmax = scale_optional(flatten_1d(raw.get("Bmax")), PHI_SIGN)
    fsa_b_over_bmin = flatten_1d(raw.get("FSA_BOverBmin"))
    fsa_b_over_bmin2 = flatten_1d(raw.get("FSA_BOverBmin2"))
    effective_passing_fraction = flatten_1d(raw.get("effectivePassingFraction"))
    xi0_trapped_boundary = flatten_1d(raw.get("xi0TrappedBoundary"))

    d_v_dr = multiply_optional(vpvol, r0)
    weight_int_vol = multiply_optional(dr, vpvol, r0)
    weight_int_area = None
    if all(value is not None for value in (dr, vpvol, gr0, r2inv, bmin)):
        with np.errstate(divide="ignore", invalid="ignore"):
            weight_int_area = dr * vpvol * gr0 * r2inv / (2.0 * np.pi * bmin)

    r_outboard = r0 + r if r0 is not None else None

    return {
        "time": time,
        "r": r,
        "dr": dr,
        "R0": r0,
        "a": a_minor,
        "R_outboard": r_outboard,
        "B0": b0,
        "GR0": gr0,
        "R2inv": r2inv,
        "Bmin": bmin,
        "Bmax": bmax,
        "FSA_B_over_Bmin": fsa_b_over_bmin,
        "FSA_B_over_Bmin2": fsa_b_over_bmin2,
        "effective_passing_fraction": effective_passing_fraction,
        "xi0_trapped_boundary": xi0_trapped_boundary,
        "phi_tor": phi_tor,
        "rho_tor": rho_tor,
        "rho_tor_norm": rho_tor_norm,
        "dVdr": d_v_dr,
        "weight_int_area": weight_int_area,
        "weight_int_vol": weight_int_vol,
        "Z0": scalar_first(raw.get("Z0_eq"), 0.0),
        "theta_eq": flatten_1d(raw.get("theta_eq")),
        "RMinusR0": raw.get("RMinusR0"),
        "ZMinusZ0": raw.get("ZMinusZ0"),
        "RMinusR0_f": raw.get("RMinusR0_f"),
        "ZMinusZ0_f": raw.get("ZMinusZ0_f"),
        "ion_Z": flatten_1d(raw.get("ion_Z")) if raw.get("ion_Z") is not None else flatten_1d(raw.get("ionmeta_Z")),
        "ion_isotopes": flatten_1d(raw.get("ion_isotopes")),
        "ion_names": parse_string_list(raw.get("ion_names") if raw.get("ion_names") is not None else raw.get("ionmeta_names")),
        "ion_hydrogennames": parse_string_list(raw.get("ion_hydrogennames")),
        "ion_tritiumnames": parse_string_list(raw.get("ion_tritiumnames")),
        "SPIMolarFraction": flatten_1d(raw.get("SPIMolarFraction")),
        "spi_Ninj": flatten_1d(raw.get("spi_Ninj")),
        "spi_xp_init": flatten_1d(raw.get("spi_xp_init")),
        "spi_vp_init": flatten_1d(raw.get("spi_vp_init")),
        "spi_rp_init": flatten_1d(raw.get("spi_rp_init")),
        "spi_ZsDrift": flatten_1d(raw.get("spi_ZsDrift")),
        "spi_isotopesDrift": flatten_1d(raw.get("spi_isotopesDrift")),
        "code": {
            "commit": raw.get("code_commit"),
            "changes": raw.get("code_changes"),
            "datetime_commit": raw.get("code_datetime_commit"),
            "datetime_simulation": raw.get("code_datetime_simulation"),
            "refspec": raw.get("code_refspec"),
        },
    }


def build_dynamic_data(raw: dict[str, Any], static: dict[str, Any]) -> dict[str, Any]:
    nt = len(static["time"])
    n_e = raw.get("n_tot") if raw.get("n_tot") is not None else raw.get("n_cold")

    dynamic = {
        "T_cold": time_radial_aligned(raw.get("T_cold"), nt),
        "n_e": time_radial_aligned(n_e, nt),
        "n_cold": time_radial_aligned(raw.get("n_cold"), nt),
        "n_i": time_aligned(raw.get("n_i"), nt),
        "n_re": time_radial_aligned(raw.get("n_re"), nt),
        "j_tot": time_radial_aligned(raw.get("j_tot"), nt),
        "j_ohm": time_radial_aligned(raw.get("j_ohm"), nt),
        "j_re": time_radial_aligned(raw.get("j_re"), nt),
        "E_field": time_radial_aligned(raw.get("E_field"), nt),
        "I_p": scale_optional(scalar_time_trace(raw.get("I_p"), nt), PHI_SIGN),
        "psi_p": scale_optional(time_radial_aligned(raw.get("psi_p"), nt), PSI_COCOS),
        "W_cold": time_radial_aligned(raw.get("W_cold"), nt),
        "W_i": time_aligned(raw.get("W_i"), nt),
        "conductivity": time_radial_aligned(raw.get("conductivity"), nt),
        "Zeff": time_radial_aligned(raw.get("Zeff"), nt),
        "runawayRate": time_radial_aligned(raw.get("runawayRate"), nt),
        "gammaDreicer": time_radial_aligned(raw.get("gammaDreicer"), nt),
        "gammaHottail": time_radial_aligned(raw.get("gammaHottail"), nt),
        "gammaTritium": time_radial_aligned(raw.get("gammaTritium"), nt),
        "gammaCompton": time_radial_aligned(raw.get("gammaCompton"), nt),
        "pCrit": scale_optional(time_radial_aligned(raw.get("pCrit"), nt), P_NORM),
        "pCritHottail": scale_optional(time_radial_aligned(raw.get("pCritHottail"), nt), P_NORM),
        "EDreic": time_radial_aligned(raw.get("EDreic"), nt),
        "Ectot": time_radial_aligned(raw.get("Ectot"), nt),
        "Tcold_radiation": time_radial_aligned(raw.get("Tcold_radiation"), nt),
        "Tcold_radiation_raw": raw.get("Tcold_radiation"),
        "x_p": reshape_spi_vector(raw.get("x_p"), nt),
        "v_p": reshape_spi_vector(raw.get("v_p"), nt),
        "Y_p": reshape_spi_scalar(raw.get("Y_p"), nt),
    }

    dynamic["I_ohm"] = current_from_j_trace(dynamic["j_ohm"], static, nt)
    dynamic["I_re"] = current_from_j_trace(dynamic["j_re"], static, nt)
    dynamic["E_cold"] = volume_integral_trace(dynamic["W_cold"], static, nt)
    dynamic["E_ion"] = volume_integral_trace(dynamic["W_i"], static, nt)
    dynamic["N_re"] = volume_integral_trace(dynamic["n_re"], static, nt)
    dynamic["n_e_volume_average"] = volume_average_trace(dynamic["n_e"], static, nt)
    dynamic["T_cold_volume_average"] = volume_average_trace(dynamic["T_cold"], static, nt)
    dynamic["Zeff_volume_average"] = volume_average_trace(dynamic["Zeff"], static, nt)

    return dynamic


def build_ids(
    raw: dict[str, Any],
    selected: list[str] | None = None,
    dd_version: str | None = None,
) -> list[Any]:
    static, dynamic = build_imas_data(raw)
    factory = make_factory(dd_version)
    ids_names = DEFAULT_IDS if selected is None else selected
    mappers = {
        "summary": map_summary,
        "runaway_electrons": map_runaway_electrons,
        "plasma_profiles": map_plasma_profiles,
        "equilibrium": map_equilibrium,
        "radiation": map_radiation,
        "spi": map_spi,
    }

    ids_list = []
    for name in ids_names:
        if name not in mappers:
            raise ValueError(f"IDS mapping is not implemented yet: {name}")
        ids = mappers[name](factory, static, dynamic)
        if ids is not None:
            ids_list.append(ids)
    return ids_list


def map_summary(factory: Any, static: dict[str, Any], dynamic: dict[str, Any]) -> Any | None:
    summary = make_ids(factory, "summary")
    if summary is None:
        return None

    time = static["time"]
    fill_ids_field(summary, "time", time)
    fill_code_metadata(summary, static)

    if static.get("R0") is not None:
        fill_ids_field(summary, "global_quantities/r0/value", static.get("R0"))
        fill_ids_field(summary, "global_quantities/r0/source", "DREAM /grid/R0")
    if static.get("B0") is not None:
        fill_ids_field(summary, "global_quantities/b0/value", np.full(time.shape, float(static["B0"])))
        fill_ids_field(summary, "global_quantities/b0/source", "derived from DREAM /grid/geometry/GR0")
    fill_ids_field(summary, "global_quantities/ip/value", dynamic.get("I_p"))
    fill_ids_field(summary, "global_quantities/current_ohm/value", dynamic.get("I_ohm"))

    e_thermal = add_optional(dynamic.get("E_cold"), dynamic.get("E_ion"))
    fill_ids_field(summary, "global_quantities/energy_electrons_thermal/value", dynamic.get("E_cold"))
    fill_ids_field(summary, "global_quantities/energy_ion_total_thermal/value", dynamic.get("E_ion"))
    fill_ids_field(summary, "global_quantities/energy_thermal/value", e_thermal)

    fill_ids_field(summary, "volume_average/t_e/value", dynamic.get("T_cold_volume_average"))
    fill_ids_field(summary, "volume_average/n_e/value", dynamic.get("n_e_volume_average"))
    fill_ids_field(summary, "volume_average/zeff/value", dynamic.get("Zeff_volume_average"))

    local_profiles = {
        "t_e": dynamic.get("T_cold"),
        "n_e": dynamic.get("n_e"),
        "zeff": dynamic.get("Zeff"),
        "e_field_parallel": dynamic.get("E_field"),
    }
    for target, data in local_profiles.items():
        fill_ids_field(summary, f"local/magnetic_axis/{target}/value", radial_sample_trace(data, 0))
        fill_ids_field(summary, f"local/separatrix/{target}/value", radial_sample_trace(data, -1))

    fill_ids_field(summary, "runaways/current/value", dynamic.get("I_re"))
    i_re = dynamic.get("I_re")
    if i_re is not None and np.asarray(i_re).size > 0:
        fill_ids_field(summary, "runaways/current_phi_max/value", float(np.nanmax(np.abs(i_re))))
    fill_ids_field(summary, "runaways/particles/value", dynamic.get("N_re"))

    return summary


def map_runaway_electrons(factory: Any, static: dict[str, Any], dynamic: dict[str, Any]) -> Any | None:
    re_ids = make_ids(factory, "runaway_electrons")
    if re_ids is None:
        return None

    time = static["time"]
    nt = len(time)
    fill_ids_field(re_ids, "time", time)
    fill_vacuum_toroidal_field(re_ids, static)

    if not resize_child_aos(re_ids, "profiles_1d", nt):
        return re_ids

    quantities = {
        "density": dynamic.get("n_re"),
        "current_density": dynamic.get("j_re"),
        "ddensity_dt_total": dynamic.get("runawayRate"),
        "ddensity_dt_dreicer": dynamic.get("gammaDreicer"),
        "ddensity_dt_hot_tail": dynamic.get("gammaHottail"),
        "ddensity_dt_tritium": dynamic.get("gammaTritium"),
        "ddensity_dt_compton": dynamic.get("gammaCompton"),
        "momentum_critical_avalanche": dynamic.get("pCrit"),
        "momentum_critical_hot_tail": dynamic.get("pCritHottail"),
        "e_field_dreicer": dynamic.get("EDreic"),
        "e_field_critical": dynamic.get("Ectot"),
    }

    for it in range(nt):
        profile = re_ids.profiles_1d[it]
        fill_1d_grid(profile, static, dynamic, it)
        for target, values in quantities.items():
            fill_ids_field(profile, target, time_slice(values, it))

    fill_ids_field(re_ids, "global_quantities/current_phi", dynamic.get("I_re"))
    return re_ids


def map_plasma_profiles(factory: Any, static: dict[str, Any], dynamic: dict[str, Any]) -> Any | None:
    pp = make_ids(factory, "plasma_profiles")
    if pp is None:
        return None

    time = static["time"]
    nt = len(time)
    fill_ids_field(pp, "time", time)
    fill_vacuum_toroidal_field(pp, static)

    if not resize_child_aos(pp, "profiles_1d", nt):
        return pp

    profiles = {
        "electrons/temperature": dynamic.get("T_cold"),
        "electrons/density": dynamic.get("n_e"),
        "electrons/density_thermal": dynamic.get("n_cold"),
        "conductivity_parallel": dynamic.get("conductivity"),
        "e_field/parallel": dynamic.get("E_field"),
        "zeff": dynamic.get("Zeff"),
        "j_total": dynamic.get("j_tot"),
        "j_ohmic": dynamic.get("j_ohm"),
    }
    ion_context = prepare_ion_profile_context(static, dynamic, nt)

    for it in range(nt):
        profile = pp.profiles_1d[it]
        fill_1d_grid(profile, static, dynamic, it)
        for target, values in profiles.items():
            fill_ids_field(profile, target, time_slice(values, it))
        fill_ion_profiles_1d(profile, static, dynamic, it, ion_context)

    return pp


def map_equilibrium(factory: Any, static: dict[str, Any], dynamic: dict[str, Any]) -> Any | None:
    eq = make_ids(factory, "equilibrium")
    if eq is None:
        return None

    time = static["time"]
    nt = len(time)
    fill_ids_field(eq, "time", time)
    fill_vacuum_toroidal_field(eq, static)

    if not resize_child_aos(eq, "time_slice", nt):
        return eq

    psi_p = dynamic.get("psi_p")
    j_tot = dynamic.get("j_tot")
    ip = dynamic.get("I_p")
    r_boundary, z_boundary = boundary_outline(static)
    profiles_2d_context = prepare_flux_surface_profiles_2d(static, psi_p)

    b_average = multiply_optional(static.get("Bmin"), static.get("FSA_B_over_Bmin"))
    gm1 = None
    if static.get("R2inv") is not None and static.get("R0") not in (None, 0.0):
        gm1 = np.asarray(static["R2inv"], dtype=float) / float(static["R0"]) ** 2
    gm5 = None
    if static.get("Bmin") is not None and static.get("FSA_B_over_Bmin2") is not None:
        gm5 = np.asarray(static["Bmin"], dtype=float) ** 2 * np.asarray(static["FSA_B_over_Bmin2"], dtype=float)
    trapped_fraction = None
    if static.get("effective_passing_fraction") is not None:
        trapped_fraction = np.clip(1.0 - np.asarray(static["effective_passing_fraction"], dtype=float), 0.0, 1.0)

    for it in range(nt):
        ts = eq.time_slice[it]
        fill_ids_field(ts, "global_quantities/ip", time_slice(ip, it))
        fill_ids_field(ts, "profiles_1d/rho_tor", static.get("rho_tor"))
        fill_ids_field(ts, "profiles_1d/rho_tor_norm", static.get("rho_tor_norm"))
        fill_ids_field(ts, "profiles_1d/phi", static.get("phi_tor"))
        fill_ids_field(ts, "profiles_1d/r_outboard", static.get("R_outboard"))
        fill_ids_field(ts, "profiles_1d/b_field_min", static.get("Bmin"))
        fill_ids_field(ts, "profiles_1d/b_field_max", static.get("Bmax"))
        fill_ids_field(ts, "profiles_1d/b_field_average", b_average)
        fill_ids_field(ts, "profiles_1d/trapped_fraction", trapped_fraction)
        fill_ids_field(ts, "profiles_1d/gm1", gm1)
        fill_ids_field(ts, "profiles_1d/gm5", gm5)
        fill_ids_field(ts, "profiles_1d/j_parallel", time_slice(j_tot, it))

        psi_arr = time_slice(psi_p, it)
        if psi_arr is not None:
            psi_arr = np.asarray(psi_arr, dtype=float)
            fill_ids_field(ts, "profiles_1d/psi", psi_arr)
            fill_ids_field(ts, "profiles_1d/psi_norm", normalized_radius(psi_arr))
            if psi_arr.size > 0:
                fill_ids_field(ts, "global_quantities/psi_magnetic_axis", float(psi_arr[0]))
                fill_ids_field(ts, "global_quantities/psi_boundary", float(psi_arr[-1]))
            d_v_dpsi = dVdpsi_from_dVdr(static.get("dVdr"), static.get("dr"), psi_arr)
            fill_ids_field(ts, "profiles_1d/dvolume_dpsi", d_v_dpsi)
            fill_equilibrium_profiles_2d(ts, psi_arr, profiles_2d_context)

        if r_boundary is not None and z_boundary is not None:
            fill_ids_field(ts, "boundary/outline/r", r_boundary)
            fill_ids_field(ts, "boundary/outline/z", z_boundary)
            fill_ids_field(ts, "boundary/geometric_axis/r", float(0.5 * (np.min(r_boundary) + np.max(r_boundary))))
            fill_ids_field(ts, "boundary/geometric_axis/z", float(0.5 * (np.min(z_boundary) + np.max(z_boundary))))

    return eq


def map_radiation(factory: Any, static: dict[str, Any], dynamic: dict[str, Any]) -> Any | None:
    radiation = make_ids(factory, "radiation")
    if radiation is None:
        return None

    emissivity_raw = dynamic.get("Tcold_radiation_raw")
    if emissivity_raw is None:
        return radiation
    emissivity = np.asarray(emissivity_raw, dtype=float)
    if emissivity.ndim != 2:
        return radiation

    time = trace_time_axis(emissivity, static["time"])
    if time is None:
        return radiation

    weights = volume_weights(static)
    if weights is None or weights.size != emissivity.shape[1]:
        return radiation

    total_power = np.sum(emissivity * weights[None, :], axis=1)
    power_inside = np.cumsum(emissivity * weights[None, :], axis=1)

    fill_ids_field(radiation, "time", time)
    if not resize_child_aos(radiation, "process", 1):
        return radiation

    process = radiation.process[0]
    fill_ids_field(process, "identifier/index", -1)
    fill_ids_field(process, "identifier/name", "dream_total_electron_radiation")
    fill_ids_field(
        process,
        "identifier/description",
        "Total cold-electron radiated power density from DREAM Tcold_radiation.",
    )

    if resize_child_aos(process, "global_quantities", time.size):
        for it in range(time.size):
            fill_ids_field(process.global_quantities[it], "inside_vessel/power_electrons", float(total_power[it]))

    if resize_child_aos(process, "profiles_1d", time.size):
        for it in range(time.size):
            profile = process.profiles_1d[it]
            fill_ids_field(profile, "grid/rho_tor_norm", static.get("rho_tor_norm"))
            fill_ids_field(profile, "grid/rho_tor", static.get("rho_tor"))
            fill_ids_field(profile, "electrons/emissivity", emissivity[it])
            fill_ids_field(profile, "electrons/power_inside", power_inside[it])

    return radiation


def map_spi(factory: Any, static: dict[str, Any], dynamic: dict[str, Any]) -> Any | None:
    spi = make_ids(factory, "spi")
    if spi is None:
        return None

    time = static["time"]
    nt = len(time)
    r0 = float(static.get("R0") or 0.0)
    fill_ids_field(spi, "time", time)

    x_p = spi_dynamic_or_init_vector(dynamic.get("x_p"), static.get("spi_xp_init"), nt)
    v_p = spi_dynamic_or_init_vector(dynamic.get("v_p"), static.get("spi_vp_init"), nt)
    y_p = spi_dynamic_or_init_scalar(dynamic.get("Y_p"), static.get("spi_rp_init"), nt)

    if x_p is None and y_p is None:
        return spi

    n_shards = x_p.shape[1] if x_p is not None else y_p.shape[1]
    if n_shards == 0:
        return spi

    groups = spi_injector_groups(static, n_shards)
    volume = None
    if y_p is not None:
        radius = np.maximum(y_p, 0.0) ** (3.0 / 5.0)
        volume = (4.0 / 3.0) * np.pi * radius**3
    components = spi_components_from_static(static)
    initial_volume = volume[0] if volume is not None and volume.shape[0] > 0 else None
    component_weights = spi_group_component_weights(groups, len(components), initial_volume)

    if not resize_child_aos(spi, "injector", len(groups)):
        return spi

    for igroup, group in enumerate(groups):
        shard_indices = group.shard_indices
        injector = spi.injector[igroup]
        fill_ids_field(injector, "name", "DREAM_SPI" if len(groups) == 1 else f"DREAM_SPI_{igroup + 1}")
        fill_ids_field(injector, "description", "Shattered pellet injection reconstructed from DREAM shard state.")
        fill_spi_core(injector, components, group, component_weights, igroup, initial_volume)

        if x_p is not None and v_p is not None and int(np.max(shard_indices)) < min(x_p.shape[1], v_p.shape[1]):
            vr0, vphi0, vz0 = dream_spi_velocity_to_imas_rphiz(x_p[0, shard_indices, :], v_p[0, shard_indices, :], r0)
            weights = shard_volume_weights(initial_volume, shard_indices)
            fill_ids_field(injector, "velocity_mass_centre_fragments_r", weighted_mean(vr0, weights))
            fill_ids_field(injector, "velocity_mass_centre_fragments_phi", weighted_mean(vphi0, weights))
            fill_ids_field(injector, "velocity_mass_centre_fragments_z", weighted_mean(vz0, weights))

        if not resize_child_aos(injector, "fragment", shard_indices.size):
            continue

        for local_index, shard_index in enumerate(shard_indices):
            fragment = injector.fragment[local_index]
            ish = int(shard_index)
            if x_p is not None and ish < x_p.shape[1]:
                r, phi, z = dream_spi_position_to_imas_rphiz(x_p[:, ish, :], r0)
                fill_ids_field(fragment, "position/r", r)
                fill_ids_field(fragment, "position/phi", phi)
                fill_ids_field(fragment, "position/z", z)
            if x_p is not None and v_p is not None and ish < x_p.shape[1] and ish < v_p.shape[1]:
                vr, vphi, vz = dream_spi_velocity_to_imas_rphiz(x_p[:, ish, :], v_p[:, ish, :], r0)
                fill_ids_field(fragment, "velocity_r", vr)
                fill_ids_field(fragment, "velocity_phi", vphi)
                fill_ids_field(fragment, "velocity_z", vz)
            elif v_p is not None and ish < v_p.shape[1]:
                fill_ids_field(fragment, "velocity_r", v_p[:, ish, 0])
                fill_ids_field(fragment, "velocity_phi", np.zeros(nt))
                fill_ids_field(fragment, "velocity_z", v_p[:, ish, 1])
            if volume is not None and ish < volume.shape[1]:
                fill_ids_field(fragment, "volume", volume[:, ish])

    return spi


def fill_vacuum_toroidal_field(ids: Any, static: dict[str, Any]) -> None:
    fill_ids_field(ids, "vacuum_toroidal_field/r0", static.get("R0"))
    b0 = static.get("B0")
    if b0 is not None:
        fill_ids_field(ids, "vacuum_toroidal_field/b0", np.full(static["time"].shape, float(b0)))


def fill_1d_grid(profile: Any, static: dict[str, Any], dynamic: dict[str, Any], it: int) -> None:
    fill_ids_field(profile, "grid/rho_tor", static.get("rho_tor"))
    fill_ids_field(profile, "grid/rho_tor_norm", static.get("rho_tor_norm"))

    psi_p = dynamic.get("psi_p")
    psi_arr = time_slice(psi_p, it)
    if psi_arr is None:
        return

    psi_arr = np.asarray(psi_arr, dtype=float)
    fill_ids_field(profile, "grid/psi", psi_arr)
    if psi_arr.size > 0:
        fill_ids_field(profile, "grid/psi_magnetic_axis", float(psi_arr[0]))
        fill_ids_field(profile, "grid/psi_boundary", float(psi_arr[-1]))
        fill_ids_field(profile, "grid/rho_pol_norm", normalized_radius(psi_arr))


def fill_code_metadata(ids: Any, static: dict[str, Any]) -> None:
    code = static.get("code", {})
    values = {key: value for key, value in code.items() if value not in (None, "")}
    if not values:
        return

    fill_ids_field(ids, "code/name", "DREAM")
    fill_ids_field(ids, "code/description", "Disruption Runaway Electron Analysis Model simulation")
    fill_ids_field(ids, "code/commit", values.get("commit"))
    fill_ids_field(ids, "code/repository", "https://github.com/chalmersplasmatheory/DREAM")
    fill_ids_field(ids, "code/version", values.get("refspec"))
    fill_ids_field(ids, "code/parameters", dream_code_parameters_xml(values))


def xml_escape(value: Any) -> str:
    text = str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def dream_code_parameters_xml(values: dict[str, Any]) -> str:
    lines = ["<dream_code_metadata>"]
    for key, value in values.items():
        if value is None:
            continue
        lines.append(f"  <{key}>{xml_escape(value)}</{key}>")
    lines.append("</dream_code_metadata>")
    return "\n".join(lines)


def cell_center_to_edges(values: Any) -> np.ndarray | None:
    arr = flatten_1d(values)
    if arr is None or arr.size < 2:
        return None
    edges = np.empty(arr.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (arr[:-1] + arr[1:])
    edges[0] = arr[0] - 0.5 * (arr[1] - arr[0])
    edges[-1] = arr[-1] + 0.5 * (arr[-1] - arr[-2])
    return edges


def dVdpsi_from_dVdr(d_v_dr: Any, dr: Any, psi: Any) -> np.ndarray | None:
    if d_v_dr is None or dr is None or psi is None:
        return None
    d_v_dr = flatten_1d(d_v_dr)
    dr = flatten_1d(dr)
    psi = flatten_1d(psi)
    psi_edges = cell_center_to_edges(psi)
    if d_v_dr is None or dr is None or psi is None or psi_edges is None:
        return None
    if d_v_dr.size != dr.size or dr.size != psi.size:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        return d_v_dr * dr / np.diff(psi_edges)


def prepare_flux_surface_profiles_2d(static: dict[str, Any], psi_p: Any) -> FluxSurface2DContext | None:
    if psi_p is None:
        return None
    psi = np.asarray(psi_p)
    if psi.ndim < 2 or psi.shape[1] == 0:
        return None

    r_minus = static.get("RMinusR0_f")
    z_minus = static.get("ZMinusZ0_f")
    if r_minus is None or z_minus is None:
        r_minus = static.get("RMinusR0")
        z_minus = static.get("ZMinusZ0")
    if r_minus is None or z_minus is None:
        return None

    r_minus = np.asarray(r_minus, dtype=float)
    z_minus = np.asarray(z_minus, dtype=float)
    if r_minus.ndim != 2 or z_minus.ndim != 2 or r_minus.shape != z_minus.shape:
        return None

    if r_minus.shape[1] == psi.shape[1] + 1:
        psi_on_edges = True
    elif r_minus.shape[1] == psi.shape[1]:
        psi_on_edges = False
    else:
        return None

    theta = static.get("theta_eq")
    if theta is None or np.asarray(theta).size != r_minus.shape[0]:
        theta = np.linspace(0.0, 2.0 * np.pi, r_minus.shape[0], endpoint=False)
    else:
        theta = np.asarray(theta, dtype=float)

    r0 = float(static.get("R0") or 0.0)
    z0 = float(static.get("Z0") or 0.0)
    return FluxSurface2DContext(
        psi_on_edges=psi_on_edges,
        theta=theta,
        r_2d=(r0 + r_minus).T,
        z_2d=(z0 + z_minus).T,
    )


def fill_equilibrium_profiles_2d(ts: Any, psi: np.ndarray, context: FluxSurface2DContext | None) -> None:
    if context is None or not resize_child_aos(ts, "profiles_2d", 1):
        return
    psi_grid = cell_center_to_edges(psi) if context.psi_on_edges else np.asarray(psi, dtype=float).reshape(-1)
    if psi_grid is None or psi_grid.size != context.r_2d.shape[0]:
        return

    p2d = ts.profiles_2d[0]
    fill_ids_field(p2d, "type/name", "total")
    fill_ids_field(p2d, "type/index", 0)
    fill_ids_field(p2d, "type/description", "Total fields")
    fill_ids_field(p2d, "grid_type/name", "inverse_psi_polar")
    fill_ids_field(p2d, "grid_type/index", 13)
    fill_ids_field(p2d, "grid_type/description", "Flux-surface grid with psi as radial label and DREAM poloidal angle as dim2.")
    fill_ids_field(p2d, "grid/dim1", psi_grid)
    fill_ids_field(p2d, "grid/dim2", context.theta)
    fill_ids_field(p2d, "r", context.r_2d)
    fill_ids_field(p2d, "z", context.z_2d)
    fill_ids_field(p2d, "psi", np.tile(psi_grid[:, np.newaxis], (1, context.theta.size)))


def boundary_outline(static: dict[str, Any]) -> tuple[np.ndarray | None, np.ndarray | None]:
    r_minus = static.get("RMinusR0_f")
    z_minus = static.get("ZMinusZ0_f")
    if r_minus is None or z_minus is None:
        r_minus = static.get("RMinusR0")
        z_minus = static.get("ZMinusZ0")
    if r_minus is None or z_minus is None:
        return None, None
    r_minus = np.asarray(r_minus, dtype=float)
    z_minus = np.asarray(z_minus, dtype=float)
    if r_minus.ndim != 2 or z_minus.ndim != 2 or r_minus.shape != z_minus.shape or r_minus.shape[1] == 0:
        return None, None
    r_outline = float(static.get("R0") or 0.0) + r_minus[:, -1]
    z_outline = float(static.get("Z0") or 0.0) + z_minus[:, -1]
    if r_outline.size > 1 and (r_outline[0] != r_outline[-1] or z_outline[0] != z_outline[-1]):
        r_outline = np.r_[r_outline, r_outline[0]]
        z_outline = np.r_[z_outline, z_outline[0]]
    return r_outline, z_outline


def reshape_spi_vector(data: Any, nt: int) -> np.ndarray | None:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 3 and arr.shape[-1] == 3:
        return arr
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    if arr.ndim == 2 and arr.shape[1] % 3 == 0:
        return arr.reshape((arr.shape[0], arr.shape[1] // 3, 3))
    return None


def reshape_spi_scalar(data: Any, nt: int) -> np.ndarray | None:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    if arr.ndim == 2:
        return arr
    return None


def spi_dynamic_or_init_vector(dynamic_value: Any, init_value: Any, nt: int) -> np.ndarray | None:
    if dynamic_value is not None:
        return dynamic_value
    init = flatten_1d(init_value)
    if init is None or init.size % 3 != 0:
        return None
    return np.repeat(init.reshape(1, init.size // 3, 3), nt, axis=0)


def spi_dynamic_or_init_scalar(dynamic_value: Any, radius_init: Any, nt: int) -> np.ndarray | None:
    if dynamic_value is not None:
        return dynamic_value
    radius = flatten_1d(radius_init)
    if radius is None:
        return None
    return np.repeat((radius ** (5.0 / 3.0)).reshape(1, -1), nt, axis=0)


def dream_spi_position_to_imas_rphiz(xyz: Any, r0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xyz = np.asarray(xyz, dtype=float)
    x = xyz[..., 0]
    z = xyz[..., 1]
    y_tor = xyz[..., 2] if xyz.shape[-1] > 2 else np.zeros_like(x)
    r_cart = r0 + x
    radius = np.hypot(r_cart, y_tor)
    phi = np.arctan2(y_tor, r_cart)
    return radius, phi, z


def dream_spi_velocity_to_imas_rphiz(xyz: Any, vxyz: Any, r0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xyz = np.asarray(xyz, dtype=float)
    vxyz = np.asarray(vxyz, dtype=float)
    x = xyz[..., 0]
    y_tor = xyz[..., 2] if xyz.shape[-1] > 2 else np.zeros_like(x)
    vx = vxyz[..., 0]
    vz = vxyz[..., 1]
    vy_tor = vxyz[..., 2] if vxyz.shape[-1] > 2 else np.zeros_like(vx)
    r_cart = r0 + x
    radius = np.hypot(r_cart, y_tor)
    radius_safe = np.where(radius > 0.0, radius, 1.0)
    velocity_r = (r_cart * vx + y_tor * vy_tor) / radius_safe
    velocity_phi = (r_cart * vy_tor - y_tor * vx) / radius_safe
    return velocity_r, velocity_phi, vz


def spi_injector_groups(static: dict[str, Any], n_shards: int) -> list[SpiInjectorGroup]:
    fractions = static.get("SPIMolarFraction")
    components = spi_components_from_static(static)
    n_components = len(components)
    if fractions is None or n_components == 0 or n_shards == 0:
        return [SpiInjectorGroup(np.arange(n_shards, dtype=int), None)]

    arr = np.asarray(fractions, dtype=float).reshape(-1)
    expected = n_components * n_shards
    if arr.size == expected + 1 and arr[0] < 0:
        arr = arr[1:]
    if arr.size != expected:
        return [SpiInjectorGroup(np.arange(n_shards, dtype=int), None)]

    shard_fractions = arr.reshape(n_components, n_shards).T
    rounded = np.round(shard_fractions, 9)
    groups: list[SpiInjectorGroup] = []
    start = 0
    for ish in range(1, n_shards):
        if not np.array_equal(rounded[ish], rounded[ish - 1]):
            groups.append(SpiInjectorGroup(np.arange(start, ish, dtype=int), np.mean(shard_fractions[start:ish], axis=0)))
            start = ish
    groups.append(SpiInjectorGroup(np.arange(start, n_shards, dtype=int), np.mean(shard_fractions[start:n_shards], axis=0)))
    valid_groups = [
        group
        for group in groups
        if group.shard_indices.size > 0 and group.fractions is not None and np.any(np.abs(group.fractions) > 1e-12)
    ]
    return valid_groups or [SpiInjectorGroup(np.arange(n_shards, dtype=int), None)]


def spi_components_from_static(static: dict[str, Any]) -> list[SpiComponent]:
    atoms = static.get("spi_Ninj")
    z_values = static.get("spi_ZsDrift")
    isotopes = static.get("spi_isotopesDrift")
    all_names = static.get("ion_names", [])
    injected_names = [(index, name) for index, name in enumerate(all_names) if "_inj" in name]
    names = [name for _, name in injected_names]
    ion_z = static.get("ion_Z")
    ion_isotopes = static.get("ion_isotopes")

    n_components = 0
    for values in (atoms, z_values, isotopes):
        if values is not None:
            n_components = max(n_components, int(np.asarray(values).size))
    if n_components == 0:
        n_components = len(names)

    components = []
    for index in range(n_components):
        z = None
        if z_values is not None and np.asarray(z_values).size > 0:
            z = int(np.asarray(z_values).reshape(-1)[min(index, np.asarray(z_values).size - 1)])
        elif index < len(injected_names) and ion_z is not None:
            name_index = injected_names[index][0]
            if name_index < np.asarray(ion_z).size:
                z = int(np.asarray(ion_z).reshape(-1)[name_index])
        if z is None:
            continue
        isotope = None
        if isotopes is not None and np.asarray(isotopes).size > 0:
            raw_isotope = int(np.asarray(isotopes).reshape(-1)[min(index, np.asarray(isotopes).size - 1)])
            isotope = raw_isotope if raw_isotope > 0 else None
        elif index < len(injected_names) and ion_isotopes is not None:
            name_index = injected_names[index][0]
            if name_index < np.asarray(ion_isotopes).size:
                raw_isotope = int(np.asarray(ion_isotopes).reshape(-1)[name_index])
                isotope = raw_isotope if raw_isotope > 0 else None
        name = names[index] if index < len(names) else element_label_from_z(z)
        atom_count = None
        if atoms is not None and index < np.asarray(atoms).size:
            atom_count = float(np.asarray(atoms).reshape(-1)[index])
        components.append(SpiComponent(name=display_species_name(name, z, isotope), z=z, isotope=isotope, atoms_n=atom_count))
    return components


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
        if initial_volume is not None and group.shard_indices.size > 0:
            max_index = int(np.max(group.shard_indices))
            if initial_volume.size > max_index:
                shard_weight = np.asarray(initial_volume[group.shard_indices], dtype=float)
        for icomp in range(n_components):
            weights[igroup, icomp] = float(np.sum(shard_weight * max(float(group.fractions[icomp]), 0.0)))
    return weights


def shard_volume_weights(initial_volume: np.ndarray | None, shard_indices: np.ndarray) -> np.ndarray | None:
    if initial_volume is None or shard_indices.size == 0:
        return None
    max_index = int(np.max(shard_indices))
    if initial_volume.size <= max_index:
        return None
    weights = np.asarray(initial_volume[shard_indices], dtype=float)
    if weights.ndim != 1 or weights.size != shard_indices.size:
        return None
    if not np.any(np.isfinite(weights)) or float(np.nansum(weights)) <= 0.0:
        return None
    return weights


def weighted_mean(values: Any, weights: np.ndarray | None) -> float:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return float("nan")
    if weights is None:
        return float(np.nanmean(values))
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return float(np.nanmean(values))
    return float(np.sum(values[valid] * weights[valid]) / np.sum(weights[valid]))


def fill_spi_core(
    injector: Any,
    components: list[SpiComponent],
    group: SpiInjectorGroup,
    component_weights: np.ndarray,
    group_index: int,
    initial_volume: np.ndarray | None,
) -> None:
    if not hasattr(injector, "pellet") or not hasattr(injector.pellet, "core"):
        return
    core = injector.pellet.core

    active: list[tuple[SpiComponent, float, float | None]] = []
    if group.fractions is None:
        for component in components:
            active.append((component, 1.0, component.atoms_n))
    else:
        component_totals = np.sum(component_weights, axis=0)
        for icomp, component in enumerate(components):
            fraction = float(group.fractions[icomp])
            if abs(fraction) <= 1e-12:
                continue
            atoms = None
            if component.atoms_n is not None and icomp < component_totals.size and component_totals[icomp] > 0.0:
                atoms = float(component.atoms_n * component_weights[group_index, icomp] / component_totals[icomp])
            active.append((component, fraction, atoms))

    atoms = [atoms_n for _, _, atoms_n in active if atoms_n is not None]
    total_atoms = float(np.sum(atoms)) if atoms and np.sum(atoms) > 0.0 else None
    if atoms:
        fill_ids_field(core, "atoms_n", total_atoms)

    total_volume = None
    if initial_volume is not None and group.shard_indices.size > 0:
        max_index = int(np.max(group.shard_indices))
        if initial_volume.size > max_index:
            volume_sum = float(np.sum(initial_volume[group.shard_indices]))
            if volume_sum > 0.0:
                total_volume = volume_sum

    if not resize_child_aos(core, "species", len(active)):
        return
    fraction_sum = None
    if group.fractions is not None:
        fraction_sum = float(np.sum([max(fraction, 0.0) for _, fraction, _ in active]))

    for index, (component, fraction, atoms_n) in enumerate(active):
        species = core.species[index]
        fill_ids_field(species, "name", component.name)
        fill_ids_field(species, "z_n", float(component.z))
        if component.isotope is not None:
            fill_ids_field(species, "a", float(component.isotope))
        density = None
        if atoms_n is not None and total_volume is not None:
            density = float(atoms_n / total_volume)
        elif total_atoms is not None and total_volume is not None and fraction_sum is not None and fraction_sum > 0.0:
            density = float(total_atoms / total_volume * max(fraction, 0.0) / fraction_sum)
        if density is not None:
            fill_ids_field(species, "density", density)


def element_label_from_z(z: int) -> str:
    labels = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 10: "Ne", 18: "Ar", 36: "Kr", 54: "Xe", 74: "W"}
    return labels.get(int(z), f"Z{int(z)}")


def display_species_name(name: str, z: int, isotope: int | None) -> str:
    if z == 1 and isotope == 1:
        return "H"
    if z == 1 and isotope == 2:
        return "D"
    if z == 1 and isotope == 3:
        return "T"
    cleaned = base_species_label(name) if "_inj" in name else name.strip()
    return cleaned or element_label_from_z(z)


def base_species_label(label: str) -> str:
    head = label.split("_", 1)[0]
    if len(head) >= 2 and head[:2].istitle():
        return head[:2]
    return head[:1]


def expected_z_from_label(label: str) -> int | None:
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
    return element_z.get(base_species_label(label))


def infer_isotope_mass_number(
    label: str,
    z: int,
    isotope: int | None,
    hydrogen_names: set[str],
    tritium_names: set[str],
) -> int | None:
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
        return 2
    if z == 10:
        return 20
    if z == 18:
        return 40
    return None


def reconcile_ion_names(
    names: list[str],
    z_list: list[int],
    isotope_list: list[int | None],
) -> list[str]:
    labels = []
    used: set[int] = set()
    for index, z in enumerate(z_list):
        chosen = None
        if index < len(names):
            z_label = expected_z_from_label(names[index])
            if z_label is None or z_label == z:
                chosen = names[index]
                used.add(index)
        if chosen is None:
            for candidate_index, name in enumerate(names):
                if candidate_index in used:
                    continue
                z_label = expected_z_from_label(name)
                if z_label is None or z_label == z:
                    chosen = name
                    used.add(candidate_index)
                    break
        if chosen is None:
            isotope = isotope_list[index] if index < len(isotope_list) else None
            chosen = f"species_{index}_Z{z}" if isotope is None else f"species_{index}_Z{z}_A{isotope}"
        labels.append(chosen)
    return labels


def prepare_ion_profile_context(
    static: dict[str, Any],
    dynamic: dict[str, Any],
    nt: int,
) -> IonProfileContext | None:
    n_i = time_aligned(dynamic.get("n_i"), nt)
    z_arr = flatten_1d(static.get("ion_Z"))
    if n_i is None or z_arr is None:
        return None
    n_i = np.asarray(n_i, dtype=float)
    if n_i.ndim != 3:
        return None

    z_list = [int(z) for z in z_arr]
    isotope_arr = flatten_1d(static.get("ion_isotopes"))
    isotope_list = []
    for index in range(len(z_list)):
        if isotope_arr is not None and index < isotope_arr.size:
            isotope = int(isotope_arr[index])
            isotope_list.append(isotope if isotope > 0 else None)
        else:
            isotope_list.append(None)

    labels = reconcile_ion_names(static.get("ion_names", []), z_list, isotope_list)
    hydrogen_names = set(static.get("ion_hydrogennames", []))
    tritium_names = set(static.get("ion_tritiumnames", []))
    a_list = [
        infer_isotope_mass_number(labels[index], z, isotope_list[index], hydrogen_names, tritium_names)
        for index, z in enumerate(z_list)
    ]
    return IonProfileContext(n_i=n_i, z_list=z_list, labels=labels, a_list=a_list)


def set_element_atomic_properties(parent: Any, z: int, mass_number: int | None) -> None:
    if not hasattr(parent, "element") or not resize_child_aos(parent, "element", 1):
        return
    element = parent.element[0]
    fill_ids_field(element, "z_n", float(z))
    if mass_number is not None:
        fill_ids_field(element, "a", float(mass_number))


def fill_ion_profiles_1d(
    profile: Any,
    static: dict[str, Any],
    dynamic: dict[str, Any],
    it: int,
    context: IonProfileContext | None,
) -> None:
    if context is None:
        return

    n_i = context.n_i
    has_ion = hasattr(profile, "ion") and resize_child_aos(profile, "ion", len(context.z_list))
    has_neutral = hasattr(profile, "neutral") and resize_child_aos(profile, "neutral", len(context.z_list))
    if not has_ion and not has_neutral:
        return

    channel_index = 0
    for species_index, z in enumerate(context.z_list):
        label = context.labels[species_index]
        mass_number = context.a_list[species_index]
        nstates = z + 1
        block = n_i[it, channel_index : min(channel_index + nstates, n_i.shape[1]), :]
        channel_index += nstates
        if block.size == 0:
            continue

        if has_ion and has_neutral:
            fill_ids_field(profile.ion[species_index], "neutral_index", species_index)
            fill_ids_field(profile.neutral[species_index], "ion_index", species_index)

        neutral_density = block[0, :]
        if has_neutral:
            neutral = profile.neutral[species_index]
            set_element_atomic_properties(neutral, z, mass_number)
            fill_ids_field(neutral, "density", neutral_density)
            if hasattr(neutral, "state") and resize_child_aos(neutral, "state", 1):
                fill_ids_field(neutral.state[0], "density", neutral_density)

        charged_block = block[1:, :]
        if has_ion:
            ion = profile.ion[species_index]
            fill_ids_field(ion, "name", label)
            set_element_atomic_properties(ion, z, mass_number)
            if charged_block.size:
                fill_ids_field(ion, "density", np.sum(charged_block, axis=0))
            if hasattr(ion, "state") and resize_child_aos(ion, "state", charged_block.shape[0]):
                for z0 in range(1, block.shape[0]):
                    state = ion.state[z0 - 1]
                    fill_ids_field(state, "name", f"{label} Z0={z0}")
                    fill_ids_field(state, "density", block[z0, :])


def parse_string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        text = value.replace("\x00", "").strip()
        if ";" in text:
            return [item.strip() for item in text.split(";") if item.strip()]
        return [item for item in text.replace(",", " ").split() if item]
    arr = np.asarray(value)
    if arr.dtype.kind in {"S", "U", "O"}:
        chars = []
        for item in arr.reshape(-1):
            text = item.decode(errors="ignore") if isinstance(item, bytes) else str(item)
            chars.append(text)
        return parse_string_list("".join(chars))
    if arr.dtype.kind in {"u", "i"}:
        text = bytes(arr.reshape(-1).astype(np.uint8)).decode(errors="ignore")
        return parse_string_list(text)
    return []


def flatten_1d(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value)
    if arr.size == 0:
        return arr.astype(float)
    return np.asarray(arr, dtype=float).reshape(-1)


def scalar_first(value: Any, default: float | None = None) -> float | None:
    if value is None:
        return default
    arr = np.asarray(value)
    if arr.size == 0:
        return default
    return float(arr.reshape(-1)[0])


def scale_optional(value: Any, factor: float) -> Any:
    if value is None:
        return None
    return np.asarray(value, dtype=float) * factor


def multiply_optional(*values: Any) -> Any:
    if any(value is None for value in values):
        return None
    result = np.asarray(values[0], dtype=float)
    for value in values[1:]:
        result = result * value
    return result


def add_optional(left: Any, right: Any) -> Any:
    if left is not None and right is not None:
        return np.asarray(left) + np.asarray(right)
    if left is not None:
        return left
    return right


def time_aligned(data: Any, nt: int) -> np.ndarray | None:
    if data is None:
        return None
    arr = np.asarray(data, dtype=float)
    if arr.ndim == 0:
        return np.full((nt,), float(arr))
    if arr.shape[0] == nt:
        return arr
    if arr.shape[0] == nt - 1:
        return np.concatenate([arr, np.take(arr, [-1], axis=0)], axis=0)
    if arr.shape[0] > nt:
        return arr[:nt]
    if 0 < arr.shape[0] < nt:
        pad_count = nt - arr.shape[0]
        pad = np.repeat(np.take(arr, [-1], axis=0), pad_count, axis=0)
        return np.concatenate([arr, pad], axis=0)
    return arr


def time_radial_aligned(data: Any, nt: int) -> np.ndarray | None:
    if data is None:
        return None
    arr = np.asarray(data, dtype=float)
    if arr.ndim == 2 and arr.shape[0] != nt and arr.shape[1] == nt:
        arr = arr.T
    return time_aligned(arr, nt)


def scalar_time_trace(data: Any, nt: int) -> np.ndarray | None:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2 and arr.shape[1] == 1:
        return arr[:, 0]
    return None


def normalized_radius(radius: Any) -> np.ndarray | None:
    if radius is None:
        return None
    arr = np.asarray(radius, dtype=float)
    if arr.size == 0:
        return arr
    denom = float(arr[-1] - arr[0])
    if abs(denom) < 1e-30:
        return np.zeros_like(arr)
    return (arr - arr[0]) / denom


def volume_weights(static: dict[str, Any]) -> np.ndarray | None:
    d_v_dr = static.get("dVdr")
    dr = static.get("dr")
    if d_v_dr is None or dr is None:
        return None
    weights = np.asarray(d_v_dr, dtype=float) * np.asarray(dr, dtype=float)
    if weights.ndim != 1 or weights.size == 0 or not np.any(np.isfinite(weights)):
        return None
    return weights


def volume_average_trace(data: Any, static: dict[str, Any], nt: int) -> np.ndarray | None:
    arr = radial_profile_trace(data, nt)
    weights = volume_weights(static)
    if arr is None or weights is None or arr.shape[1] != weights.size:
        return None
    denom = np.sum(weights)
    if denom == 0:
        return None
    return np.sum(arr * weights[None, :], axis=1) / denom


def volume_integral_trace(data: Any, static: dict[str, Any], nt: int) -> np.ndarray | None:
    aligned = time_aligned(data, nt)
    weights = volume_weights(static)
    if aligned is None or weights is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == weights.size:
        return np.sum(arr * weights[None, :], axis=1)
    if arr.ndim == 3 and arr.shape[-1] == weights.size:
        return np.sum(arr * weights[None, None, :], axis=(1, 2))
    return None


def trace_time_axis(data: np.ndarray, full_time: np.ndarray) -> np.ndarray | None:
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


def current_from_j_trace(data: Any, static: dict[str, Any], nt: int) -> np.ndarray | None:
    arr = radial_profile_trace(data, nt)
    weights = static.get("weight_int_area")
    if arr is None or weights is None:
        return None
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or weights.size != arr.shape[1]:
        return None
    return np.sum(arr * weights[None, :], axis=1)


def radial_profile_trace(data: Any, nt: int) -> np.ndarray | None:
    aligned = time_aligned(data, nt)
    if aligned is None:
        return None
    arr = np.asarray(aligned, dtype=float)
    if arr.ndim == 2:
        return arr
    return None


def radial_sample_trace(data: Any, index: int) -> np.ndarray | None:
    if data is None:
        return None
    arr = np.asarray(data, dtype=float)
    if arr.ndim != 2 or arr.shape[1] == 0:
        return None
    return arr[:, index]


def time_slice(data: Any, index: int) -> Any:
    if data is None:
        return None
    arr = np.asarray(data)
    if arr.ndim == 0:
        return arr.item()
    if index >= arr.shape[0]:
        return None
    return arr[index]
