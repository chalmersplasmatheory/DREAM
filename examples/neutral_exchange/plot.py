#!/usr/bin/env python3
"""Plot neutral temperatures and check conservation in both species orders."""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import e

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-show', action='store_true', help='Save without opening a window.')
    args = parser.parse_args()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    curves = {}

    for label, order, style in (
        ('forward', ('D', 'Ne'), '-'),
        ('reverse', ('Ne', 'D'), '--'),
    ):
        with h5py.File(HERE / f'output_{label}.h5', 'r') as f:
            t = f['grid/t'][:].reshape(-1)
            W = f['eqsys/W_n'][:].reshape(t.size, 2)
            density = f['eqsys/n_i'][:].reshape(t.size, 13)

        # n_i contains each species' charge states consecutively.
        neutral_columns = [0, 2 if order[0] == 'D' else 11]
        n = density[:, neutral_columns]
        T = W / (1.5 * e * n)
        total = W.sum(axis=1)
        relative_error = (total - total[0]) / total[0]
        equilibrium = np.sum(n[0] * T[0]) / n[0].sum()

        assert np.all(np.isfinite(T)), f'{label}: non-finite temperature'
        assert np.all(W >= 0), f'{label}: negative energy'
        np.testing.assert_allclose(n, np.broadcast_to(n[0], n.shape), rtol=1e-12)
        assert np.max(np.abs(relative_error)) < 1e-8, f'{label}: energy not conserved'
        np.testing.assert_allclose(T[-1], equilibrium, rtol=1e-3)
        assert np.all(np.diff(T[:, order.index('D')]) >= -1e-8)
        assert np.all(np.diff(T[:, order.index('Ne')]) <= 1e-8)

        curves[label] = T[:, order.index('D')]
        for name, color in (('D', 'tab:blue'), ('Ne', 'tab:orange')):
            axes[0].plot(t * 1e3, T[:, order.index(name)], style,
                         color=color, label=f'{name}, {label}')
        axes[1].plot(t * 1e3, relative_error, style, label=label)
        print(f'{label}: final temperatures {dict(zip(order, T[-1]))} eV; '
              f'max relative energy error = {np.max(np.abs(relative_error)):.3e}')

    difference = np.max(np.abs(curves['forward'] - curves['reverse']))
    print(f'Max D temperature difference after swapping species order: {difference:.6g} eV')
    print('Species-order invariance: ' + ('PASS' if difference < 1e-6 else 'FAIL (known coefficient issue)'))
    axes[0].axhline(6, color='gray', linewidth=1, label='Expected equilibrium: 6 eV')
    axes[0].set(xlabel='Time [ms]', ylabel='Neutral temperature [eV]',
                title='Neutral D–Ne equilibration', xlim=(0, 20))
    axes[1].set(xlabel='Time [ms]', ylabel='(Total energy − initial) / initial',
                title='Neutral energy conservation')
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.savefig(HERE / 'neutral_exchange.png', dpi=180)
    fig.savefig(HERE / 'neutral_exchange.pdf')
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
