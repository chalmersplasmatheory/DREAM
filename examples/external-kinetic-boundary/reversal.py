#!/usr/bin/env python3
"""Demonstrate two-way exchange at the runaway grid's outer momentum face."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.constants import c, e

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'py'))

from DREAM import DREAMOutput, DREAMSettings
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.DistributionFunction as Dist
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


P_MIN = 1.0
P_MAX = 12.0
N_P = 88
N_XI = 28
SEED = 1e15  # m^-3
FIELD_TIMES = np.array([0.0, 2.1e-3, 2.3e-3, 5.8e-3])
FIELD_VALUES = np.array([1.0, 1.0, -1.0, -1.0])  # V/m
NSTEPS = 116


def settings(enabled):
    ds = DREAMSettings()
    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setPmin(P_MIN)
    ds.runawaygrid.setPmax(P_MAX)
    ds.runawaygrid.setNp(N_P)
    ds.runawaygrid.setNxi(N_XI)

    ds.radialgrid.setB0(5)
    ds.radialgrid.setMinorRadius(0.22)
    ds.radialgrid.setWallRadius(0.22)
    ds.radialgrid.setNr(1)
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_ULTRA_RELATIVISTIC
    ds.eqsys.n_i.addIon('D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=1e19)
    ds.eqsys.T_cold.setPrescribedData(1000)
    ds.eqsys.E_field.setPrescribedData(FIELD_VALUES, times=FIELD_TIMES)

    # Start a passing packet near p_max. The cylindrical momentum measure
    # is 2*pi*p^2 dp dxi. DREAM initializes n_re from the complete f_re.
    p_faces = np.linspace(P_MIN, P_MAX, N_P+1)
    xi_faces = np.linspace(-1, 1, N_XI+1)
    p = (p_faces[:-1] + p_faces[1:]) / 2
    xi = (xi_faces[:-1] + xi_faces[1:]) / 2
    f = np.exp(-0.5*((p[None, :] - 10.95)/0.45)**2 - 0.5*((xi[:, None] - 0.86)/0.075)**2)
    seed_density = 2*np.pi * np.sum(f * (p**2 * np.diff(p_faces))[None, :] * np.diff(xi_faces)[:, None])
    dist = ds.eqsys.f_re
    dist.setInitialValue(f[None, :, :] * (SEED / seed_density), r=[0], p=p, xi=xi)
    dist.setAdvectionInterpolationMethod(ad_int=Dist.AD_INTERP_UPWIND)
    dist.setBoundaryCondition(Dist.BC_F_0)
    dist.setExternalKineticBoundary(enabled)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.timestep.setTmax(FIELD_TIMES[-1])
    ds.timestep.setNt(NSTEPS)
    ds.output.setFilename('output.h5')
    if enabled:
        ds.other.include('fluid/gammaExternal', 'fluid/gammaExternalIn', 'fluid/gammaExternalOut',
                         'fluid/gammaExternalPos', 'fluid/gammaExternalNeg')
    return ds


def run_case(binary, case, enabled):
    case.mkdir(parents=True)
    settings(enabled).save(str(case / 'settings.h5'))
    with (case / 'run.log').open('w') as log:
        result = subprocess.run([str(binary), 'settings.h5'], cwd=case, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"DREAM exited with code {result.returncode}. See {case / 'run.log'}")


def read_case(case, enabled):
    output = DREAMOutput(str(case / 'output.h5'), lazy=False, loadsettings=False)
    eq = output.eqsys
    t = output.grid.t
    if (t.size != NSTEPS+1 or np.any(np.diff(t) <= 0)
            or not np.isclose(t[-1], FIELD_TIMES[-1], rtol=1e-12, atol=0)):
        raise ValueError(f"{case}: the time advance is incomplete.")
    field = eq.E_field[:, 0]
    np.testing.assert_allclose(field, np.interp(t, FIELD_TIMES, FIELD_VALUES), rtol=0, atol=1e-12)
    result = dict(t=t, field=field, total=eq.n_re[:, 0], resolved=eq.f_re.density()[:, 0],
                  f=eq.f_re[:, 0, :, :], p=eq.f_re.momentumgrid.p1, xi=eq.f_re.momentumgrid.p2)
    if enabled:
        result['positive'] = eq.n_external_pos[:, 0]
        result['negative'] = eq.n_external_neg[:, 0]
        result['current'] = eq.j_re[:, 0]
        result['resolved_current'] = eq.f_re.currentDensity()[:, 0]
        diagnostics = output.other.fluid
        for name in ['gammaExternal', 'gammaExternalIn', 'gammaExternalOut', 'gammaExternalPos', 'gammaExternalNeg']:
            result[name] = diagnostics[name][:, 0]
    output.close()
    if not all(np.isfinite(values).all() for values in result.values()):
        raise ValueError(f"{case}: nonfinite output.")
    return result


def check(two_way, legacy):
    t, p, n = (two_way[key] for key in ['t', 'positive', 'negative'])
    dt = np.diff(t)
    # Other quantities are sampled at the accepted step ends in DREAM.
    outward = two_way['gammaExternalOut'] * dt
    inward = two_way['gammaExternalIn'] * dt
    cumulative_out = np.r_[0, np.cumsum(outward)]
    cumulative_in = np.r_[0, np.cumsum(inward)]
    transfer_error = np.max(np.abs(np.vstack([
        np.diff(p) - dt*two_way['gammaExternalPos'],
        np.diff(n) - dt*two_way['gammaExternalNeg'],
        (p+n)[1:] - (p[0]+n[0]) - cumulative_out[1:] + cumulative_in[1:],
        dt*two_way['gammaExternal'] - outward + inward,
    ]))) / SEED
    balance_error = np.max(np.abs(two_way['resolved'] + p + n - two_way['total'])) / SEED
    drift = np.max(np.abs(two_way['total'] - two_way['total'][0])) / SEED
    current_error = np.max(np.abs(two_way['current'] - two_way['resolved_current'] - e*c*(p-n))) / (e*c*SEED)
    np.testing.assert_allclose(two_way['total'][0], legacy['total'][0], rtol=1e-12)
    np.testing.assert_allclose(two_way['f'][0], legacy['f'][0], rtol=1e-12)
    if max(transfer_error, balance_error, drift, current_error) > 1e-7 or min(p.min(), n.min()) < -1e-10*SEED:
        raise ValueError("The two-way run failed its signed particle or current balance.")
    fill = t[1:] <= FIELD_TIMES[1]
    reversed_field = t[1:] >= FIELD_TIMES[2]
    if outward[fill].sum() < 0.1*SEED or inward[reversed_field].sum() < 0.01*SEED:
        raise ValueError("The case must show appreciable outward transfer and return after field reversal.")
    if two_way['resolved'][-1] <= legacy['resolved'][-1]:
        raise ValueError("The two-way run should retain more resolved electrons than the legacy control.")

    metrics = dict(
        initial_density_m3=float(two_way['total'][0]),
        max_external_fraction=float(np.max(p+n)/SEED),
        outward_before_ramp_fraction=float(outward[fill].sum()/SEED),
        inward_after_ramp_fraction=float(inward[reversed_field].sum()/SEED),
        two_way_final_resolved_fraction=float(two_way['resolved'][-1]/SEED),
        legacy_final_resolved_fraction=float(legacy['resolved'][-1]/SEED),
        max_normalized_density_budget_error=float(balance_error),
        max_normalized_flux_budget_error=float(transfer_error),
        max_normalized_total_drift=float(drift),
        max_normalized_current_error=float(current_error),
    )
    return metrics, cumulative_out, cumulative_in


def save_figure(fig, directory, name):
    fig.savefig(directory / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(directory / f'{name}.png', bbox_inches='tight', dpi=200)


def plot(two_way, legacy, cumulative_out, cumulative_in, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.rcParams.update({
        'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
        'legend.frameon': False, 'pdf.fonttype': 42,
    })
    t = two_way['t']*1e3
    blue, orange, grey, green = '#0072B2', '#C45214', '#525D66', '#007B64'
    p, n = two_way['positive']/SEED, two_way['negative']/SEED
    ramp = FIELD_TIMES[1:3]*1e3

    fig, (field, budget, transfer) = plt.subplots(3, 1, figsize=(9.2, 7.0), sharex=True,
        layout='constrained', gridspec_kw={'height_ratios': [1, 2.1, 1.8]})
    fig.suptitle('Density budget during an electric-field reversal', fontsize=16, fontweight='bold')
    for ax in (field, budget, transfer):
        ax.axvspan(*ramp, color='#737F8D', alpha=0.14, lw=0)
        ax.grid(axis='y', alpha=0.2)
        ax.set_xlim(0, t[-1])
    field.plot(t, two_way['field'], color=grey, lw=2)
    field.axhline(0, color=grey, alpha=0.4, lw=1)
    field.set(ylabel=r'$E$ (V/m)', yticks=[-1, 0, 1])
    field.annotate('Field reverses', (np.mean(ramp), 0), xytext=(ramp[1]+0.35, 0.32),
                   arrowprops={'arrowstyle': '-', 'color': grey}, color=grey)

    budget.plot(t, two_way['total']/SEED, color=grey, lw=1.5, ls=':', label='Total runaway density')
    budget.plot(t, two_way['resolved']/SEED, color=blue, lw=2.4, label='Resolved, two-way boundary')
    budget.plot(t, legacy['resolved']/SEED, color=blue, lw=1.5, ls='--', label='Resolved, legacy boundary')
    budget.plot(t, p, color=orange, lw=2, label=r'External, $\xi>0$')
    budget.plot(t, n, color=green, lw=1.6, label=r'External, $\xi<0$')
    budget.set(ylabel='Density / initial density', ylim=(-0.025, 1.08))
    budget.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=9)

    transfer.plot(t, cumulative_out/SEED, color=orange, lw=2, label='Out of grid')
    transfer.plot(t, cumulative_in/SEED, color=blue, lw=2, label='Back into grid')
    transfer.plot(t, p+n, color=grey, lw=1.5, ls='--', label='Net external density')
    transfer.set(xlabel='Time (ms)', ylabel='Cumulative transfer\n/ initial density')
    transfer.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=9)
    save_figure(fig, directory, 'density-budget')
    plt.close(fig)

    peak_index = np.argmax(p+n)
    indices = (0, peak_index, len(t)-1)
    names = ('Initial packet', 'Peak external density', 'After reversal')
    fig, axes = plt.subplots(2, 3, figsize=(10.4, 5.8), sharex=True, sharey=True, layout='constrained')
    fig.suptitle('Resolved runaway distribution at the outer momentum boundary', fontsize=15, fontweight='bold')
    reference = two_way['f'][0].max()
    levels = np.geomspace(1e-5, 1, 21)
    for col, (i, name) in enumerate(zip(indices, names)):
        axes[0, col].set_title(f'{name}\n{t[i]:.2f} ms', fontsize=11)
        for row, result in enumerate((two_way, legacy)):
            ax = axes[row, col]
            contour = ax.contourf(result['p'], result['xi'],
                np.maximum(result['f'][i]/reference, levels[0]), levels=levels,
                norm=LogNorm(levels[0], levels[-1]), cmap='viridis', extend='max')
            ax.axvline(P_MAX, color='#182C35', lw=1.2)
            ax.set(xlim=(7, P_MAX+0.06), ylim=(0.05, 1), xticks=[7, 8, 9, 10, 11, 12])
            if col == 2:
                ax.annotate(r'$p_\mathrm{max}$', (P_MAX, 0.16), xytext=(-6, 0),
                            textcoords='offset points', ha='right', color='white', fontsize=9)
            if row == 1:
                ax.set_xlabel(r'$p/(m_e c)$')
    axes[0, 0].set_ylabel('Two-way boundary\nPitch $\\xi$')
    axes[1, 0].set_ylabel('Legacy boundary\nPitch $\\xi$')
    fig.colorbar(contour, ax=axes, label=r'$f_\mathrm{re}/\max[f_\mathrm{re}(0)]$',
                 ticks=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    save_figure(fig, directory, 'distribution-contours')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, help='New directory for the runs and figures')
    parser.add_argument('--binary', type=Path, default=REPO / 'build/iface/dreami',
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
        run_case(binary, directory / 'two-way', enabled=True)
        run_case(binary, directory / 'legacy', enabled=False)
    two_way = read_case(directory / 'two-way', enabled=True)
    legacy = read_case(directory / 'legacy', enabled=False)
    metrics, cumulative_out, cumulative_in = check(two_way, legacy)
    plot(two_way, legacy, cumulative_out, cumulative_in, directory)
    summary = json.dumps(metrics, indent=2)
    (directory / 'summary.json').write_text(summary + '\n')
    print(summary)
    print(f'Figures and logs are in {directory}')


if __name__ == '__main__':
    main()
