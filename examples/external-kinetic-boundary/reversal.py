#!/usr/bin/env python3
"""Compare both poles of the two-way kinetic boundary with a wider grid."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.constants import c, e, m_e

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'py'))

from DREAM import DREAMOutput, DREAMSettings
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.DistributionFunction as Dist
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


CUT = 1
WIDE = 2
P_CELLS_PER_UNIT = 80
N_XI = 64
NSTEPS = 400
SEED = 1e9  # m^-3
FIELD = 1.0  # V/m
# The signed momentum excursion must cross p_max=1 at both poles.
LEG = 1.06*m_e*c/(e*FIELD)
RAMP_HALF_WIDTH = 2*(4*LEG/NSTEPS)
FIELD_TIMES = np.array([
    0, LEG-RAMP_HALF_WIDTH, LEG+RAMP_HALF_WIDTH,
    3*LEG-RAMP_HALF_WIDTH, 3*LEG+RAMP_HALF_WIDTH, 4*LEG,
])
FIELD_VALUES = FIELD*np.array([1, 1, -1, -1, 1, 1])
DIAGNOSTICS = (
    'gammaExternal', 'gammaExternalOut', 'gammaExternalIn',
    'gammaExternalPos', 'gammaExternalNeg',
)


def settings(pmax, external):
    ds = DREAMSettings()
    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setPmin(0)
    ds.runawaygrid.setPmax(pmax)
    ds.runawaygrid.setNp(pmax*P_CELLS_PER_UNIT)
    ds.runawaygrid.setNxi(N_XI)

    ds.radialgrid.setB0(5)
    ds.radialgrid.setMinorRadius(0.22)
    ds.radialgrid.setWallRadius(0.22)
    ds.radialgrid.setNr(1)
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_ULTRA_RELATIVISTIC
    ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e14)
    ds.eqsys.T_cold.setPrescribedData(1000)
    ds.eqsys.E_field.setPrescribedData(FIELD_VALUES, times=FIELD_TIMES)

    # Isotropic seed centered on the origin. Both grids share the same cells
    # below p=1, so they also share the same discretized initial distribution.
    dp = 1/P_CELLS_PER_UNIT
    dxi = 2/N_XI
    p = dp*(np.arange(pmax*P_CELLS_PER_UNIT)+0.5)
    xi = -1+dxi*(np.arange(N_XI)+0.5)
    profile = np.exp(-0.5*(p/0.10)**2)
    profile[p >= 0.6] = 0
    f = np.tile(profile, (N_XI, 1))
    f *= SEED/(2*np.pi*dp*dxi*np.sum(p**2*f))
    dist = ds.eqsys.f_re
    dist.setInitialValue(f[None, :, :], r=[0], p=p, xi=xi)
    dist.setAdvectionInterpolationMethod(ad_int=Dist.AD_INTERP_UPWIND)
    dist.setBoundaryCondition(Dist.BC_F_0)
    dist.setExternalKineticBoundary(external)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.timestep.setTmax(FIELD_TIMES[-1])
    ds.timestep.setNt(NSTEPS)
    ds.output.setFilename('output.h5')
    if external:
        ds.other.include([f'fluid/{name}' for name in DIAGNOSTICS])
    return ds


def run_case(binary, case, pmax, external):
    case.mkdir(parents=True)
    settings(pmax, external).save(str(case/'settings.h5'))
    with (case/'run.log').open('w') as log:
        result = subprocess.run([str(binary), 'settings.h5'], cwd=case,
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"DREAM exited with code {result.returncode}. See {case/'run.log'}.")


def read_case(case, external):
    output = DREAMOutput(str(case/'output.h5'), lazy=False, loadsettings=False)
    eq = output.eqsys
    t = output.grid.t
    if (t.size != NSTEPS+1 or np.any(np.diff(t) <= 0)
            or not np.isclose(t[-1], FIELD_TIMES[-1], rtol=1e-12, atol=0)):
        raise ValueError(f"{case}: the time advance is incomplete.")
    field = eq.E_field[:, 0]
    np.testing.assert_allclose(field, np.interp(t, FIELD_TIMES, FIELD_VALUES), rtol=0, atol=1e-12)
    f = eq.f_re[:, 0, :, :]
    grid = eq.f_re.momentumgrid
    p, xi = grid.p1, grid.p2
    # Use the same momentum-cell measure as the DREAM density moment.
    weight = grid.Vprime_VpVol[0] * grid.DP1[0] * grid.DP2[0]/SEED
    particles = f*weight
    resolved = particles.sum(axis=(1, 2))
    result = dict(t=t, field=field, f=f, p=p, xi=xi, p_faces=grid.p1_f, resolved=resolved,
                  mean_parallel=np.sum(particles*xi[:, None]*p[None, :], axis=(1, 2))/resolved)
    if external:
        result['positive'] = eq.n_external_pos[:, 0]/SEED
        result['negative'] = eq.n_external_neg[:, 0]/SEED
        result['total'] = eq.n_re[:, 0]/SEED
        result['external_current'] = (eq.j_re[:, 0]-eq.f_re.currentDensity()[:, 0])/(e*c*SEED)
        for name in DIAGNOSTICS:
            result[name] = output.other.fluid[name][:, 0]/SEED
    else:
        beyond = particles[:, :, p > CUT]
        result['beyond_positive'] = beyond[:, xi > 0, :].sum(axis=(1, 2))
        result['beyond_negative'] = beyond[:, xi < 0, :].sum(axis=(1, 2))
    output.close()
    if not all(np.isfinite(values).all() for values in result.values()):
        raise ValueError(f"{case}: nonfinite output.")
    return result


def check(cut, wide):
    t = cut['t']
    np.testing.assert_allclose(t, wide['t'], rtol=0, atol=1e-12)
    np.testing.assert_allclose(cut['f'][0], wide['f'][0][:, wide['p'] < CUT], rtol=1e-12)
    np.testing.assert_allclose(wide['resolved'][0], 1, rtol=1e-12)
    dt = np.diff(t)
    positive, negative = cut['positive'], cut['negative']
    outward = dt*cut['gammaExternalOut']
    inward = dt*cut['gammaExternalIn']
    budget_error = np.max(np.abs(cut['resolved']+positive+negative-cut['total']))
    flux_error = np.max(np.abs(np.stack([
        np.diff(positive)-dt*cut['gammaExternalPos'],
        np.diff(negative)-dt*cut['gammaExternalNeg'],
        np.diff(positive+negative)-outward+inward,
        dt*cut['gammaExternal']-outward+inward,
    ])))
    current_error = np.max(np.abs(cut['external_current']-positive+negative))
    total_drift = np.max(np.abs(cut['total']-cut['total'][0]))
    wide_drift = np.max(np.abs(wide['resolved']-wide['resolved'][0]))
    if (max(budget_error, flux_error, current_error, total_drift) > 1e-7
            or min(positive.min(), negative.min()) < -1e-10):
        raise ValueError("The cut-grid run failed its signed particle or current balance.")
    if wide_drift > 0.001:
        raise ValueError("The wide reference lost appreciable density at p=2.")
    peaks = [int(np.argmax(positive)), int(np.argmax(negative))]
    if not t[peaks[0]] < 2*LEG < t[peaks[1]]:
        raise ValueError("The two pole peaks must occur on opposite legs of the pulse.")
    middle = int(np.argmin(np.abs(t-2*LEG)))
    stop = len(t)-1
    if (min(positive[peaks[0]], negative[peaks[1]]) < 0.2
            or positive[middle] > 0.1 or (positive+negative)[stop] > 0.12):
        raise ValueError("Both poles must load and then return most of their density.")
    if (min(wide['beyond_positive'].max(), wide['beyond_negative'].max()) < 0.2
            or np.max(np.abs(wide['mean_parallel'][[middle, stop]])) > 0.15
            or np.max(np.abs(cut['mean_parallel'][[middle, stop]])) > 0.15):
        raise ValueError("The packet must cross both cut faces and return near zero parallel momentum.")

    return dict(
        initial_density_m3=SEED,
        peak_external_positive=float(positive[peaks[0]]),
        peak_external_negative=float(negative[peaks[1]]),
        peak_wide_beyond_cut_positive=float(wide['beyond_positive'].max()),
        peak_wide_beyond_cut_negative=float(wide['beyond_negative'].max()),
        external_positive_at_midpoint=float(positive[middle]),
        final_external_positive=float(positive[stop]),
        final_external_negative=float(negative[stop]),
        final_cut_resolved=float(cut['resolved'][stop]),
        final_wide_resolved=float(wide['resolved'][stop]),
        final_cut_mean_p_parallel=float(cut['mean_parallel'][stop]),
        final_wide_mean_p_parallel=float(wide['mean_parallel'][stop]),
        max_cut_density_budget_error=float(budget_error),
        max_cut_flux_budget_error=float(flux_error),
        max_cut_total_drift=float(total_drift),
        max_cut_current_error=float(current_error),
        max_wide_resolved_drift=float(wide_drift),
    ), [0, peaks[0], middle, peaks[1], stop]


def save_figure(fig, directory, name):
    fig.savefig(directory/f'{name}.pdf', bbox_inches='tight')
    fig.savefig(directory/f'{name}.png', bbox_inches='tight', dpi=200)


def plot_budget(cut, wide, directory, plt):
    t = cut['t']*1e3
    blue, orange, grey, purple = '#0072B2', '#C45114', '#515B65', '#65478B'
    fig, (field, budget, poles) = plt.subplots(
        3, 1, figsize=(9, 6.8), sharex=True, layout='constrained',
        gridspec_kw={'height_ratios': [1, 1.65, 1.65]},
    )
    fig.suptitle('Particle exchange at both pitch poles', fontsize=16, fontweight='bold')
    for ax in (field, budget, poles):
        for i in (1, 3):
            ax.axvspan(1e3*FIELD_TIMES[i], 1e3*FIELD_TIMES[i+1], color=grey, alpha=0.13, lw=0)
        ax.grid(axis='y', alpha=0.2)
        ax.set_xlim(0, t[-1])
    field.plot(t, cut['field'], color=grey, lw=2)
    field.axhline(0, color=grey, alpha=0.4, lw=1)
    field.set(ylabel=r'$E$ (V/m)', yticks=[-1, 0, 1])
    for x, label in ((LEG, 'First reversal'), (3*LEG, 'Second reversal')):
        field.annotate(label, (1e3*x, 0), xytext=(1e3*x+0.15, 0.35), color=grey, fontsize=9)

    budget.plot(t, cut['total'], ls=':', lw=2.5, color=purple, label='Cut grid + both external poles')
    budget.plot(t, wide['resolved'], ls='--', lw=1.5, color=grey, label='Resolved, wide grid')
    budget.plot(t, cut['resolved'], lw=2.2, color=purple, label='Resolved, cut grid')
    budget.set(ylabel='Density / initial density', ylim=(0, 1.07))
    budget.legend(loc='lower center', ncols=1, fontsize=9)

    poles.plot(t, wide['beyond_positive'], ls='--', lw=1.6, color=orange, alpha=0.75,
               label=r'Wide grid, $p>1$ and $\xi>0$')
    poles.plot(t, wide['beyond_negative'], ls='--', lw=1.6, color=blue, alpha=0.75,
               label=r'Wide grid, $p>1$ and $\xi<0$')
    poles.plot(t, cut['positive'], lw=2.3, color=orange, label=r'Cut grid, external $\xi>0$')
    poles.plot(t, cut['negative'], lw=2.3, color=blue, label=r'Cut grid, external $\xi<0$')
    poles.set(xlabel='Time (ms)', ylabel='Density / initial density', ylim=(-0.018, 0.9))
    handles, labels = poles.get_legend_handles_labels()
    order = [2, 3, 0, 1]
    poles.legend([handles[i] for i in order], [labels[i] for i in order],
                 loc='upper center', ncols=2, fontsize=9)
    save_figure(fig, directory, 'density-budget')
    plt.close(fig)


def plot_distribution(cut, wide, snapshots, directory, plt):
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(2, len(snapshots), figsize=(13.2, 5.0), sharex=True,
                             sharey=True, layout='constrained')
    fig.suptitle('Two reversals in momentum space', fontsize=16, fontweight='bold')
    titles = ('Initial packet', 'Positive excursion', 'Crossing near zero', 'Negative excursion', 'End of pulse')
    reference = wide['f'][0].max()
    levels = np.geomspace(1e-4, 1, 21)
    angle = np.linspace(0, np.pi, 301)

    for row, result in enumerate((wide, cut)):
        # Extend the plotted cell averages to the domain faces. In particular,
        # omitting xi=+/-1 leaves a spurious blank wedge along p_parallel.
        p = np.r_[result['p_faces'][0], result['p'], result['p_faces'][-1]][None, :]
        xi = np.r_[-1, result['xi'], 1][:, None]
        parallel = p*xi
        perpendicular = p*np.sqrt(1-xi**2)
        for col, (i, title) in enumerate(zip(snapshots, titles)):
            ax = axes[row, col]
            if row == 0:
                ax.set_title(f"{title}\n{1e3*result['t'][i]:.2f} ms", fontsize=10)
            else:
                ax.set_facecolor('#EDF0F3')
                ax.set_xlabel(r'$p_\parallel/(m_e c)$')
            values = np.pad(result['f'][i]/reference, 1, mode='edge')
            contour = ax.contourf(parallel, perpendicular,
                np.maximum(values, levels[0]), levels=levels,
                norm=LogNorm(levels[0], levels[-1]), cmap='viridis', extend='max')
            ax.plot(CUT*np.cos(angle), CUT*np.sin(angle), color='#E3B66B', ls='--', lw=1.2)
            ax.set(xlim=(-1.65, 1.65), ylim=(0, 0.65), xticks=[-1, 0, 1], yticks=[0, 0.3, 0.6])
    axes[0, 0].set_ylabel('Wide grid, $p_\\mathrm{max}=2$\n' + r'$p_\perp/(m_e c)$')
    axes[1, 0].set_ylabel('Two-way, $p_\\mathrm{max}=1$\n' + r'$p_\perp/(m_e c)$')
    fig.colorbar(contour, ax=axes, label=r'$f_\mathrm{re}/\max[f_\mathrm{re}(0)]$',
                 ticks=[1e-4, 1e-3, 1e-2, 1e-1, 1])
    save_figure(fig, directory, 'distribution-contours')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, help='New directory for the runs and figures')
    parser.add_argument('--binary', type=Path, default=REPO/'build/iface/dreami',
                        help='Path to the compiled dreami executable')
    parser.add_argument('--analyze-only', action='store_true', help='Regenerate figures from existing output')
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not args.analyze_only:
        if directory.exists():
            parser.error(f'{directory} exists. Choose a new directory or use --analyze-only.')
        binary = args.binary.resolve()
        if not binary.is_file():
            parser.error(f'{binary} is missing. Build dreami or pass --binary.')
        run_case(binary, directory/'cut', CUT, external=True)
        run_case(binary, directory/'wide', WIDE, external=False)

    cut = read_case(directory/'cut', external=True)
    wide = read_case(directory/'wide', external=False)
    metrics, snapshots = check(cut, wide)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
        'legend.frameon': False, 'pdf.fonttype': 42,
    })
    plot_budget(cut, wide, directory, plt)
    plot_distribution(cut, wide, snapshots, directory, plt)
    summary = json.dumps(metrics, indent=2)
    (directory/'summary.json').write_text(summary + '\n')
    print(summary)
    print(f'Figures and logs are in {directory}')


if __name__ == '__main__':
    main()
