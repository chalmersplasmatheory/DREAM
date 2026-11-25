#!/usr/bin/env python3
"""Validate the isotropic finite-detrapping Rechester-Rosenbluth coefficient.

For each selected magnetic perturbation this example runs two hot-tail cases:

    1. ordinary isotropic Rechester-Rosenbluth transport,
    2. the eigenmode finite-detrapping correction.

The coefficient check uses the DREAM output coefficients `hottail/Drr`,
`hottail/Drr_detrapping` and `hottail/nu_D_f2`. The raw DREAM
detrapping-limited coefficient is compared to an independent reconstruction
from the analytic geometry and the local deflection frequency.
"""

import argparse
import logging
import os
from pathlib import Path
import sys
import warnings

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

MPLCONFIGDIR = Path(os.environ.get("TMPDIR", "/tmp")) / "dream-rr-detrapping-matplotlib"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

warnings.filterwarnings("ignore", message="'created' timestamp seems very low.*")
warnings.filterwarnings("ignore", message="'modified' timestamp seems very low.*")
logging.getLogger("fontTools.ttLib.tables._h_e_a_d").setLevel(logging.ERROR)

sys.path.insert(0, str(REPO_ROOT / "py"))

import numpy as np
from scipy.constants import c
from scipy.linalg import eigh_tridiagonal
from scipy.special import jn_zeros

from DREAM import DREAMOutput, DREAMSettings, runiface
import DREAM.Settings.AdvectionInterpolation as AdInt
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectrons as NCold
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport


MU0 = 4e-7 * np.pi
J01_SQUARED = jn_zeros(0, 1)[0]**2
MODE_RR = "rr"
MODE_EIGENMODE = "eigenmode"
MODELS = (MODE_RR, MODE_EIGENMODE)


def configure_matplotlib(usetex=True):
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.4,
        "legend.handlelength": 2.3,
        "axes.formatter.use_mathtext": True,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "figure.constrained_layout.use": True,
        "mathtext.fontset": "cm",
        "text.usetex": bool(usetex),
        "text.latex.preamble": r"\usepackage{amsmath}",
    })


def parse_dbb(values):
    """Accept either '--dBB 1e-3 1e-2' or '--dBB 1e-3,1e-2'."""
    out = []
    for value in values:
        for part in str(value).replace(",", " ").split():
            out.append(float(part))
    return out


def shape_profiles(r, args):
    """Analytic shaped-flux-surface profiles used by DREAM's radial grid."""
    r = np.asarray(r, dtype=float)
    rho2 = (r / args.minor_radius) ** 2

    return {
        "GOverR0": args.B0 * np.ones_like(r),
        "psi": -MU0 * args.plasma_current * args.minor_radius * (1.0 - rho2),
        "psi_prime": 2.0 * MU0 * args.plasma_current * r / args.minor_radius,
        "kappa": 1.0 + (args.elongation - 1.0) * rho2,
        "kappa_prime": 2.0 * (args.elongation - 1.0) * r / args.minor_radius**2,
        "delta": args.triangularity * rho2,
        "delta_prime": 2.0 * args.triangularity * r / args.minor_radius**2,
        "Delta": args.shafranov_shift * rho2,
        "Delta_prime": 2.0 * args.shafranov_shift * r / args.minor_radius**2,
    }


def configure_radial_grid(ds, args):
    a = args.minor_radius
    r_ref = np.linspace(0.0, a, 64)
    profiles = shape_profiles(r_ref, args)

    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
    ds.radialgrid.setShaping(
        psi=profiles["psi"],
        rpsi=r_ref,
        GOverR0=profiles["GOverR0"],
        rG=r_ref,
        kappa=profiles["kappa"],
        rkappa=r_ref,
        delta=profiles["delta"],
        rdelta=r_ref,
        Delta=profiles["Delta"],
        rDelta=r_ref,
    )
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a)
    ds.radialgrid.setMajorRadius(args.major_radius)
    ds.radialgrid.setNr(args.nr)


def initial_profiles(args):
    r = np.linspace(0.0, args.minor_radius, 8)
    t = np.array([0.0])

    n_hot = args.hot_fraction * args.density
    n_hot *= 1.0 - args.hollowing * (r / args.minor_radius) ** 2
    n_hot = np.maximum(n_hot, 0.05 * args.hot_fraction * args.density)

    return {
        "r": r,
        "t": t,
        "E": args.electric_field * np.ones((1, r.size)),
        "T_cold": args.cold_temperature * np.ones((1, r.size)),
        "n_hot": n_hot,
        "T_hot": args.hot_temperature * np.ones_like(r),
    }


def configure_collisions(ds):
    ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED
    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.collisions.pstar_mode = Collisions.PSTAR_MODE_COLLISIONAL
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL


def configure_plasma(ds, prof, args):
    r = prof["r"]

    ds.eqsys.E_field.setPrescribedData(prof["E"], radius=r, times=prof["t"])
    ds.eqsys.T_cold.setPrescribedData(prof["T_cold"], radius=r, times=prof["t"])
    ds.eqsys.n_cold.setType(NCold.TYPE_SELFCONSISTENT)

    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)
    ds.eqsys.n_i.setAdvectionInterpolationMethodCharged(
        ad_int=AdInt.AD_INTERP_UPWIND,
        ad_jac=AdInt.AD_INTERP_JACOBIAN_UPWIND,
        fluxlimiterdamping=1.0,
    )
    ds.eqsys.n_i.addIon(
        name="D",
        Z=1,
        iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED,
        n=args.density + prof["n_hot"],
        r=r,
    )
    if args.neon_density > 0:
        ds.eqsys.n_i.addIon(
            name="Ne",
            Z=10,
            iontype=Ions.IONS_DYNAMIC_NEUTRAL,
            n=args.neon_density * np.ones_like(r),
            r=r,
        )


def disable_runaway_sources(ds):
    ds.eqsys.n_re.setDreicer(False)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
    ds.eqsys.n_re.setCompton(False)
    ds.eqsys.n_re.setTritium(False)
    ds.eqsys.n_re.setHottail(False)
    ds.runawaygrid.setEnabled(False)


def configure_hot_tail(ds, prof, args):
    ds.hottailgrid.setEnabled(True)
    ds.hottailgrid.setNxi(1)
    ds.hottailgrid.setNp(args.np)
    ds.hottailgrid.setPmax(args.pmax)

    ds.eqsys.f_hot.setInitialProfiles(
        n0=prof["n_hot"],
        T0=prof["T_hot"],
        rn0=prof["r"],
        rT0=prof["r"],
    )
    ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0)
    ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_ZERO)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(
        ad_int=FHot.AD_INTERP_TCDF,
        ad_jac=FHot.AD_INTERP_JACOBIAN_FULL,
    )
    ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)


def configure_transport(ds, model, dBB, args):
    r = np.linspace(0.0, args.minor_radius, args.nr + 1)
    t = np.array([0.0, args.tmax])
    dBB_t_r = dBB * np.ones((t.size, r.size))

    if model == MODE_RR:
        ds.eqsys.f_hot.transport.setMagneticPerturbation(dBB_t_r, t=t, r=r)
    elif model == MODE_EIGENMODE:
        ds.eqsys.f_hot.transport.setDetrappingRechesterRosenbluth(
            dBB_t_r,
            t=t,
            r=r,
        )
    else:
        raise ValueError(f"Unknown transport model '{model}'.")


def configure_solver_and_output(ds, output_path, args):
    ds.solver.setType(Solver.LINEAR_IMPLICIT)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.solver.preconditioner.setEnabled(False)
    ds.solver.setVerbose(False)

    ds.timestep.setTmax(args.tmax)
    ds.timestep.setNt(args.nt)
    ds.output.setFilename(str(output_path))
    ds.output.setTiming(stdout=False, file=False)
    ds.other.include(
        "fluid",
        "scalar",
        "hottail/Drr",
        "hottail/Drr_detrapping",
        "hottail/kappa_detrapping",
        "hottail/nu_D_f2",
    )


def make_settings(model, dBB, args, output_path):
    ds = DREAMSettings()
    profiles = initial_profiles(args)

    configure_collisions(ds)
    configure_radial_grid(ds, args)
    configure_plasma(ds, profiles, args)
    disable_runaway_sources(ds)
    configure_hot_tail(ds, profiles, args)
    configure_transport(ds, model, dBB, args)
    configure_solver_and_output(ds, output_path, args)

    return ds


def case_record(dBB, model, outdir):
    tag = f"dBB_{dBB:.4e}_{model}"
    return {
        "dBB": dBB,
        "model": model,
        "settings": outdir / f"{tag}_settings.h5",
        "output": outdir / f"{tag}_output.h5",
    }


def run_case(dBB, model, args, outdir):
    record = case_record(dBB, model, outdir)
    ds = make_settings(model, dBB, args, record["output"])
    ds.save(str(record["settings"]))

    label = f"{model:11s} dB/B={dBB:.4e}"
    if args.no_run:
        print(f"[{label}] settings written: {record['settings']}")
    else:
        print(f"[{label}] running")
        runiface(ds, str(record["output"]), quiet=args.quiet)
        print(f"[{label}] output written: {record['output']}")

    return record


def load_results(records):
    return {
        (record["dBB"], record["model"]): DREAMOutput(str(record["output"]))
        for record in records
    }


def extract_drr(output, time_index):
    """Return `hottail/Drr` as an array with shape `(nr+1, np)`."""
    drr = np.asarray(output.other.hottail.Drr[time_index], dtype=float)
    if drr.ndim == 3:
        if drr.shape[1] != 1:
            raise RuntimeError(f"Expected isotropic Drr with one pitch cell, got {drr.shape}.")
        return drr[:, 0, :]
    if drr.ndim == 2:
        return drr
    raise RuntimeError(f"Unexpected Drr shape: {drr.shape}.")


def extract_dream_detrap(output, time_index):
    """Return DREAM's raw `hottail/Drr_detrapping` diagnostic."""
    if not hasattr(output.other.hottail, "Drr_detrapping"):
        raise RuntimeError(
            "This output does not contain hottail/Drr_detrapping. "
            "Rebuild DREAM and rerun this script."
        )

    drr = np.asarray(output.other.hottail.Drr_detrapping[time_index], dtype=float)
    if drr.ndim == 3:
        if drr.shape[1] != 1:
            raise RuntimeError(f"Expected isotropic Drr_detrapping with one pitch cell, got {drr.shape}.")
        return drr[:, 0, :]
    if drr.ndim == 2:
        return drr
    raise RuntimeError(f"Unexpected Drr_detrapping shape: {drr.shape}.")


def extract_nuD_f2(output, time_index):
    """Return `hottail/nu_D_f2` as an array with shape `(nr, np)`."""
    if not hasattr(output.other.hottail, "nu_D_f2"):
        raise RuntimeError(
            "This output does not contain hottail/nu_D_f2. "
            "Rerun this script instead of using older coefficient-only output."
        )

    nuD = np.asarray(output.other.hottail.nu_D_f2[time_index], dtype=float)
    if nuD.ndim == 3:
        # TrappingLimitedRRTransport uses the first P2 slice in isotropic PXi.
        return nuD[:, 0, :]
    if nuD.ndim == 2:
        return nuD
    raise RuntimeError(f"Unexpected nu_D_f2 shape: {nuD.shape}.")


def interpolate_to_radial_faces(cell_quantity, output):
    """Match the radial-face interpolation used by FVM DiffusionTerm."""
    r = np.asarray(output.grid.r, dtype=float)
    r_f = np.asarray(output.grid.r_f, dtype=float)
    nr = r.size

    face_quantity = np.empty((nr + 1, cell_quantity.shape[1]))
    face_quantity[0, :] = cell_quantity[0, :]
    face_quantity[nr, :] = cell_quantity[nr - 1, :]

    dr = np.diff(r_f)
    dr_f = np.diff(r)
    for ir in range(1, nr):
        w_right = 0.5 * dr[ir - 1] / dr_f[ir - 1]
        face_quantity[ir, :] = (
            w_right * cell_quantity[ir, :]
            + (1.0 - w_right) * cell_quantity[ir - 1, :]
        )

    return face_quantity


def split_legendre_quadrature(n):
    """Legendre points split at pi so Bmin/Bmax are resolved cleanly."""
    x, w = np.polynomial.legendre.leggauss(n)
    theta = []
    weight = []
    for lo, hi in [(0.0, np.pi), (np.pi, 2.0 * np.pi)]:
        theta.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo))
        weight.append(0.5 * (hi - lo) * w)
    return np.concatenate(theta), np.concatenate(weight)


def trapped_pitch_eigenvalue(xiT):
    """Lowest Dirichlet eigenvalue from a centered finite-difference pitch grid."""
    pitch = np.linspace(-xiT, xiT, 257)
    spacing = pitch[1] - pitch[0]
    faces = (pitch[1:] + pitch[:-1]) / 2
    diffusion = (1 - faces**2) / (2 * spacing**2)
    return eigh_tridiagonal(
        diffusion[:-1] + diffusion[1:], -diffusion[1:-1],
        eigvals_only=True, select='i', select_range=(0, 0)
    )[0]


def transport_geometry_factors(r_f, args, nquad=800):
    """Return the RR geometry factor, trapped boundary and eigenmode cap per unit nu_D."""
    theta, weight = split_legendre_quadrature(nquad)
    profiles = shape_profiles(r_f, args)
    geometric = np.empty_like(r_f, dtype=float)
    xiT = np.zeros_like(r_f, dtype=float)
    detrap_per_nuD = np.zeros_like(r_f, dtype=float)

    ct = np.cos(theta)
    st = np.sin(theta)

    for k, r in enumerate(r_f):
        if r <= 0:
            geometric[k] = 0.5
            continue

        delta = profiles["delta"][k]
        kappa = profiles["kappa"][k]
        Delta = profiles["Delta"][k]
        delta_prime = profiles["delta_prime"][k]
        kappa_prime = profiles["kappa_prime"][k]
        Delta_prime = profiles["Delta_prime"][k]

        cdt = np.cos(delta * st)
        sdt = np.sin(delta * st)
        sin_theta_delta = st * cdt + sdt * ct
        cos_theta_delta = ct * cdt - st * sdt
        r_kappa_prime = r * kappa_prime

        J_over_Rr = (
            kappa * (cdt + Delta_prime * ct)
            + st
            * sin_theta_delta
            * (r_kappa_prime + ct * (delta * (kappa + r_kappa_prime) - r * kappa * delta_prime))
        )
        R_over_R0 = 1.0 + (Delta + r * cos_theta_delta) / args.major_radius
        jacobian = r * R_over_R0 * J_over_Rr

        delta_term = (1.0 + delta * ct) * sin_theta_delta
        kappa_term = kappa * ct
        nabla_r2 = (kappa_term**2 + delta_term**2) / J_over_Rr**2

        B_tor = profiles["GOverR0"][k] / R_over_R0
        B_pol = profiles["psi_prime"][k] / (2.0 * np.pi * R_over_R0)
        B = np.sqrt(B_tor**2 + nabla_r2 * B_pol**2)

        Bmin = np.min(B)
        Bmax = np.max(B)
        B_over_Bmin = B / Bmin
        fsa_B = np.sum(weight * jacobian * B_over_Bmin) / np.sum(weight * jacobian)
        xiT2 = max(0.0, 1.0 - Bmin / Bmax)
        geometric[k] = 0.5 * fsa_B * (1.0 - xiT2)
        xiT[k] = np.sqrt(xiT2)

        if xiT[k] > 0:
            kappa_pitch = trapped_pitch_eigenvalue(xiT[k])
            detrap_per_nuD[k] = args.minor_radius**2 * kappa_pitch / J01_SQUARED

    return geometric, xiT, detrap_per_nuD


def rr_geometry_factor(r_f, args):
    """Mirror the analytic toroidal geometry used by DREAM for ordinary RR."""
    return transport_geometry_factors(r_f, args)[0]


def analytic_rr(output, dBB, args):
    p = np.asarray(output.grid.hottail.p, dtype=float)
    r_f = np.asarray(output.grid.r_f, dtype=float)
    v = c * p / np.sqrt(1.0 + p**2)
    geometric = rr_geometry_factor(r_f, args)
    return np.pi * args.major_radius * dBB**2 * geometric[:, None] * v[None, :]


def analytic_detrap(output, args, time_index):
    """Reconstruct the implemented eigenmode detrapping limit directly."""
    r_f = np.asarray(output.grid.r_f, dtype=float)
    nuD_cell = extract_nuD_f2(output, time_index)
    nuD_face = interpolate_to_radial_faces(nuD_cell, output)
    _, _, detrap_per_nuD = transport_geometry_factors(r_f, args)
    return detrap_per_nuD[:, None] * np.maximum(nuD_face, 0.0)


def harmonic_limit(D1, D2):
    D1, D2 = np.broadcast_arrays(D1, D2)
    out = np.zeros_like(D1, dtype=float)
    good = (D1 > 0) & (D2 > 0)
    out[good] = D1[good] * D2[good] / (D1[good] + D2[good])
    return out


def expected_eigenmode_rr(Drr, Ddetrap, trapped):
    """Apply the cap on trapped faces and ordinary RR at untrapped faces."""
    return np.where(trapped[:, None], harmonic_limit(Drr, Ddetrap), Drr)


def relative_error(value, reference, floor=1e-300):
    if np.size(value) == 0:
        return 0.0

    scale = np.maximum(np.abs(reference), floor)
    err = np.abs(value - reference) / scale
    return np.max(err) if np.all(np.isfinite(err)) else np.inf


def extract_dream_detrap_limits(results, dBB_values, time_index):
    """Read the raw D_detrap component calculated inside DREAM."""
    return {
        dBB: extract_dream_detrap(results[(dBB, MODE_EIGENMODE)], time_index)
        for dBB in dBB_values
    }


def reconstruct_detrap_limits(results, dBB_values, args):
    return {
        dBB: analytic_detrap(
            results[(dBB, MODE_EIGENMODE)],
            args,
            args.validation_index,
        )
        for dBB in dBB_values
    }


def validate_coefficients(results, dBB_values, detrap_analytic, detrap_dream, args):
    rows = []

    for dBB in dBB_values:
        rr_output = results[(dBB, MODE_RR)]
        Drr = extract_drr(rr_output, args.validation_index)
        Deff = extract_drr(results[(dBB, MODE_EIGENMODE)], args.validation_index)
        expected_rr = analytic_rr(rr_output, dBB, args)
        expected_detrap = detrap_analytic[dBB]
        dream_detrap = detrap_dream[dBB]
        eigen = results[(dBB, MODE_EIGENMODE)]
        trapped = eigen.other.hottail.kappa_detrapping[args.validation_index] > 0
        expected_eff = expected_eigenmode_rr(Drr, dream_detrap, trapped)
        rr_error = relative_error(Drr, expected_rr)
        limiter_error = relative_error(Deff[trapped], expected_eff[trapped])
        detrap_error = relative_error(dream_detrap, expected_detrap)
        axis_error = relative_error(Deff[~trapped], Drr[~trapped])

        rows.append((dBB, rr_error, limiter_error, detrap_error, axis_error))

    return rows


def hot_inventory(output):
    return np.asarray(output.eqsys.n_hot.integral(), dtype=float)


def fit_loss_time(t, inventory):
    good = (inventory > 0) & np.isfinite(inventory)
    t = t[good]
    inventory = inventory[good]
    if t.size < 4:
        return np.nan

    lo = max(1, int(0.15 * t.size))
    hi = max(lo + 3, int(0.75 * t.size))
    slope, _ = np.polyfit(t[lo:hi], np.log(inventory[lo:hi] / inventory[0]), 1)
    return np.nan if slope >= 0 else -1.0 / slope


def plot_edge_coefficients(results, dBB_values, detrap_analytic, detrap_dream, args, figdir):
    import matplotlib.pyplot as plt

    dBB = dBB_values[len(dBB_values) // 2]
    rr_output = results[(dBB, MODE_RR)]
    p = np.asarray(rr_output.grid.hottail.p, dtype=float)
    ir = -1

    Drr = extract_drr(rr_output, args.validation_index)
    Deff = extract_drr(results[(dBB, MODE_EIGENMODE)], args.validation_index)
    expected_rr = analytic_rr(rr_output, dBB, args)
    expected_detrap = detrap_analytic[dBB]
    dream_detrap = detrap_dream[dBB]
    eigen = results[(dBB, MODE_EIGENMODE)]
    trapped = eigen.other.hottail.kappa_detrapping[args.validation_index] > 0
    expected_eff = expected_eigenmode_rr(Drr, dream_detrap, trapped)
    detrap_mask = np.isfinite(expected_detrap[ir, :]) & (expected_detrap[ir, :] > 0)
    dream_mask = np.isfinite(dream_detrap[ir, :]) & (dream_detrap[ir, :] > 0)

    colors = {
        "rr": "#5b2a86",
        "eff": "#0f766e",
        "detrap": "#8c1d18",
    }

    fig, ax = plt.subplots(figsize=(5.1, 3.0))
    ax.loglog(p, Drr[ir, :], "-", color=colors["rr"], label=r"DREAM $D_{\rm RR}$")
    ax.loglog(p, expected_rr[ir, :], "--", color='black', label=r"analytic $D_{\rm RR}$")
    ax.loglog(p, Deff[ir, :], "-", color=colors["eff"], label=r"DREAM $D_{\rm eff}$")
    ax.loglog(p, expected_eff[ir, :], "--", color=colors["eff"], label=r"harmonic $D_{\rm eff}$")
    ax.loglog(
        p[detrap_mask],
        expected_detrap[ir, detrap_mask],
        ":",
        color=colors["detrap"],
        lw=1.8,
        label=r"analytic $D_{\rm detrap}$",
    )
    ax.loglog(
        p[dream_mask],
        dream_detrap[ir, dream_mask],
        ".",
        color=colors["detrap"],
        ms=2.0,
        alpha=0.55,
        label=r"DREAM $D_{\rm detrap}$",
    )
    ax.set_xlabel(r"$p/(m_ec)$")
    ax.set_ylabel(r"$D_{rr}\;[\mathrm{m^2\,s^{-1}}]$")
    ax.set_title(rf"$\delta B/B={dBB:.1e}$, outer radial face")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))

    path = figdir / "rr_eigenmode_coefficient_check.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_limiter_identity(results, dBB_values, detrap_dream, args, figdir):
    import matplotlib.pyplot as plt

    x_values = []
    rr_ratios = []
    detrap_ratios = []

    for dBB in dBB_values:
        Drr = extract_drr(results[(dBB, MODE_RR)], args.validation_index)
        Deff = extract_drr(results[(dBB, MODE_EIGENMODE)], args.validation_index)
        Ddetrap = detrap_dream[dBB]
        mask = (Drr > 0) & (Ddetrap > 0) & np.isfinite(Deff)

        x = Drr[mask] / Ddetrap[mask]
        x_values.append(x)
        rr_ratios.append(Deff[mask] / Drr[mask])
        detrap_ratios.append(Deff[mask] / Ddetrap[mask])

    x = np.concatenate(x_values)
    rr_ratios = np.concatenate(rr_ratios)
    detrap_ratios = np.concatenate(detrap_ratios)

    order = np.argsort(x)
    xline = np.logspace(np.log10(np.nanmin(x)), np.log10(np.nanmax(x)), 400)

    fig, ax = plt.subplots(figsize=(5.1, 3.3))
    ax.loglog(x[order], rr_ratios[order], ".", ms=2.0, alpha=0.45, label=r"$D_{\rm eff}/D_{\rm RR}$")
    ax.loglog(x[order], detrap_ratios[order], ".", ms=2.0, alpha=0.45, label=r"$D_{\rm eff}/D_{\rm detrap}$")
    ax.loglog(xline, 1.0 / (1.0 + xline), "k-", lw=1.0)
    ax.loglog(xline, xline / (1.0 + xline), "k--", lw=1.0)
    ax.axvline(1.0, color="0.5", ls=":", lw=1.0)
    ax.set_xlabel(r"$D_{\rm RR}/D_{\rm detrap}$")
    ax.set_ylabel("limiter ratio")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)

    path = figdir / "rr_eigenmode_limiter_identity.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_rr_vs_detrap(results, dBB_values, detrap_dream, args, figdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.8, 3.35))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(dBB_values)))

    finite_values = []
    for color, dBB in zip(colors, dBB_values):
        Drr = extract_drr(results[(dBB, MODE_RR)], args.validation_index)
        Ddetrap = detrap_dream[dBB]
        mask = (Drr > 0) & (Ddetrap > 0) & np.isfinite(Ddetrap)
        ax.loglog(
            Drr[mask],
            Ddetrap[mask],
            ".",
            ms=2.5,
            alpha=0.45,
            color=color,
            label=rf"$\delta B/B={dBB:.1e}$",
        )
        finite_values.append(Drr[mask])
        finite_values.append(Ddetrap[mask])

    values = np.concatenate(finite_values)
    lo = np.nanmin(values)
    hi = np.nanmax(values)
    line = np.logspace(np.log10(lo), np.log10(hi), 300)
    ax.loglog(line, line, "k:", lw=1.0, label=r"$D_{\rm RR}=D_{\rm detrap}$")

    ax.set_xlabel(r"$D_{\rm RR}\;[\mathrm{m^2\,s^{-1}}]$")
    ax.set_ylabel(r"$D_{\rm detrap}\;[\mathrm{m^2\,s^{-1}}]$")
    ax.legend(frameon=False, fontsize=7)

    path = figdir / "rr_vs_detrap_coefficient.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_detrap_validation(results, dBB_values, detrap_analytic, detrap_dream, args, figdir):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.8, 3.35))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(dBB_values)))

    finite_values = []
    for color, dBB in zip(colors, dBB_values):
        Ddetrap = detrap_analytic[dBB]
        dream = detrap_dream[dBB]
        mask = (Ddetrap > 0) & np.isfinite(Ddetrap) & np.isfinite(dream)
        if not np.any(mask):
            continue

        ax.loglog(
            Ddetrap[mask],
            dream[mask],
            ".",
            ms=2.5,
            alpha=0.45,
            color=color,
            label=rf"$\delta B/B={dBB:.1e}$",
        )
        finite_values.append(Ddetrap[mask])
        finite_values.append(dream[mask])

    values = np.concatenate(finite_values)
    lo = np.nanmin(values)
    hi = np.nanmax(values)
    line = np.logspace(np.log10(lo), np.log10(hi), 300)
    ax.loglog(line, line, "k:", lw=1.0, label=r"exact agreement")

    ax.set_xlabel(r"analytic $D_{\rm detrap}\;[\mathrm{m^2\,s^{-1}}]$")
    ax.set_ylabel(r"DREAM $D_{\rm detrap}\;[\mathrm{m^2\,s^{-1}}]$")
    ax.legend(frameon=False, fontsize=7)

    path = figdir / "rr_eigenmode_detrap_validation.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_hot_inventory(results, dBB_values, args, figdir):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(4.8, 3.35))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(dBB_values)))

    for color, dBB in zip(colors, dBB_values):
        for model, linestyle in [(MODE_RR, "-"), (MODE_EIGENMODE, "--")]:
            output = results[(dBB, model)]
            t = np.asarray(output.grid.t, dtype=float)
            N = hot_inventory(output)
            ax.plot(1e6 * t, N / N[0], color=color, ls=linestyle)

    ax.set_xlabel(r"$t\;[\mu\mathrm{s}]$")
    ax.set_ylabel(r"$N_{\rm hot}(t)/N_{\rm hot}(0)$")
    ax.set_ylim(bottom=0)

    handles = [
        Line2D([0], [0], color="0.1", ls="-", marker="o", label="RR"),
        Line2D([0], [0], color="0.1", ls="--", marker="s", label="eigenmode"),
    ]
    handles += [
        Line2D([0], [0], color=color, lw=1.8, label=rf"{dBB:.1e}")
        for color, dBB in zip(colors, dBB_values)
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=min(6, len(handles)),
        fontsize=7,
    )

    path = figdir / "rr_eigenmode_nhot_inventory.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def write_summary(rows, figure_paths, outdir):
    lines = [
        "# RR Eigenmode Detrapping Validation",
        "",
        "| dB/B | RR max rel. err. | limiter max rel. err. | detrap max rel. err. | axis rel. err. |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for dBB, rr_error, limiter_error, detrap_error, axis_error in rows:
        lines.append(
            f"| {dBB:.6e} | {rr_error:.6e} | {limiter_error:.6e} | "
            f"{detrap_error:.6e} | {axis_error:.6e} |"
        )

    lines += ["", "Figures:"]
    lines += [f"- {path.name}" for path in figure_paths if path is not None]

    path = outdir / "validation_summary.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ordinary and finite-detrapping RR hot-tail transport checks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--dBB", nargs="+", default=["3e-4", "3e-3", "3e-2"],
        help="Magnetic perturbation values. Comma-separated lists are accepted.")
    parser.add_argument("--output-dir", default=str(SCRIPT_DIR / "output"),
        help="Directory for settings, outputs and figures.")

    parser.add_argument("--nr", type=int, default=10, help="Number of radial cells.")
    parser.add_argument("--np", type=int, default=90, help="Number of momentum cells.")
    parser.add_argument("--pmax", type=float, default=2.5, help="Maximum hot-tail momentum.")
    parser.add_argument("--tmax", type=float, default=3.0e-4, help="Simulation time [s].")
    parser.add_argument("--nt", type=int, default=80, help="Number of time steps.")

    parser.add_argument("--minor-radius", type=float, default=0.5, help="Minor radius [m].")
    parser.add_argument("--major-radius", type=float, default=3.0, help="Major radius [m].")
    parser.add_argument("--B0", type=float, default=5.0, help="On-axis toroidal magnetic field [T].")
    parser.add_argument("--plasma-current", type=float, default=8.0e5,
        help="Reference plasma current for the analytic poloidal flux [A].")
    parser.add_argument("--elongation", type=float, default=1.55,
        help="Edge elongation of the analytic toroidal flux surfaces.")
    parser.add_argument("--triangularity", type=float, default=0.0,
        help="Edge triangularity of the analytic toroidal flux surfaces.")
    parser.add_argument("--shafranov-shift", type=float, default=0.0,
        help="Edge Shafranov shift of the analytic toroidal flux surfaces [m].")

    parser.add_argument("--density", type=float, default=5.0e19, help="Cold electron density [m^-3].")
    parser.add_argument("--neon-density", type=float, default=1.0e19,
        help="Initial flat neutral neon density [m^-3]. Set to 0 to disable neon.")
    parser.add_argument("--cold-temperature", type=float, default=1.0e3, help="Cold electron temperature [eV].")
    parser.add_argument("--hot-temperature", type=float, default=7.5e4, help="Initial f_hot temperature [eV].")
    parser.add_argument("--hot-fraction", type=float, default=5.0e-3, help="Initial hot density fraction.")
    parser.add_argument("--hollowing", type=float, default=0.65, help="Quadratic radial hollowing of n_hot.")
    parser.add_argument("--electric-field", type=float, default=0.0, help="Prescribed electric field [V/m].")

    parser.add_argument("--tolerance", type=float, default=5.0e-3,
        help="Relative tolerance shown in figures and used for validation.")
    parser.add_argument("--validation-index", type=int, default=0,
        help="Drr time index used for coefficient validation.")
    parser.add_argument("--no-run", action="store_true", help="Only write settings files.")
    parser.add_argument("--plot-only", action="store_true", help="Use existing outputs and regenerate plots.")
    parser.add_argument("--skip-plots", action="store_true", help="Run cases but do not generate figures.")
    parser.add_argument("--quiet", action="store_true", help="Run DREAM quietly.")
    parser.add_argument("--no-usetex", action="store_false", dest="usetex",
        help="Disable external LaTeX text rendering.")
    parser.set_defaults(usetex=True)
    return parser.parse_args()


def build_records(dBB_values, outdir):
    return [case_record(dBB, model, outdir) for dBB in dBB_values for model in MODELS]


def main():
    args = parse_args()
    dBB_values = parse_dbb(args.dBB)
    outdir = Path(args.output_dir).resolve()
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    records = build_records(dBB_values, outdir)
    if not args.plot_only:
        records = [run_case(dBB, model, args, outdir) for dBB in dBB_values for model in MODELS]

    if args.no_run:
        return 0

    missing = [record["output"] for record in records if not record["output"].exists()]
    if missing:
        raise SystemExit("Missing output files:\n" + "\n".join(str(path) for path in missing))

    results = load_results(records)
    detrap_dream = extract_dream_detrap_limits(
        results,
        dBB_values,
        time_index=args.validation_index,
    )
    detrap_analytic = reconstruct_detrap_limits(results, dBB_values, args)
    summary = validate_coefficients(
        results,
        dBB_values,
        detrap_analytic,
        detrap_dream,
        args,
    )

    figures = []
    if not args.skip_plots:
        configure_matplotlib(usetex=args.usetex)
        figures = [
            plot_edge_coefficients(results, dBB_values, detrap_analytic, detrap_dream, args, figdir),
            plot_limiter_identity(results, dBB_values, detrap_dream, args, figdir),
            plot_rr_vs_detrap(results, dBB_values, detrap_dream, args, figdir),
            plot_detrap_validation(results, dBB_values, detrap_analytic, detrap_dream, args, figdir),
            plot_hot_inventory(results, dBB_values, args, figdir),
        ]
    summary_path = write_summary(summary, figures, outdir)

    print("\nValidation summary")
    for dBB, rr_error, limiter_error, detrap_error, axis_error in summary:
        print(
            f"  dB/B={dBB:.4e}: RR err={rr_error:.3e}, "
            f"limiter err={limiter_error:.3e}, "
            f"detrap err={detrap_error:.3e}, axis err={axis_error:.3e}"
        )
    print(f"  figures: {figdir}")
    print(f"  summary: {summary_path}")

    failed = [row for row in summary if max(row[1:]) > args.tolerance]
    if failed:
        raise SystemExit(
            "Coefficient validation exceeded tolerance. "
            f"See {summary_path} and {figdir}."
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
