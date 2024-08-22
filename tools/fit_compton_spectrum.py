#!/usr/bin/env python3
#

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit
import sys


def func(E, lG, c1, c2, c3):
    """
    Function to fit to photon spectrum.
    """
    z = (np.log(E) + c1) / c2 + c3 * E**2
    return lG - np.exp(-z)-z+1


def eval_integral(lG, c1, c2, c3):
    """
    Evaluate the integral of the photon spectrum from 0 to inf.
    """
    I, _ = quad(lambda E : np.exp(func(E, lG, c1, c2, c3)), 0, np.inf)
    return I


def fit_spectrum(E, G):
    """
    Fit curve to the provided photon spectrum.
    """
    popt, _ = curve_fit(func, E, np.log(G), p0=(np.amax(G), 1.2, 0.8, 0), bounds=([0, -np.inf, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf]))
    return popt


def load_csv(filename): return load_delim(filename, ',')
def load_tsv(filename): return load_delim(filename, '\t')


def load_delim(filename, delim):
    """
    Load photon spectrum from a delimited text file.
    """
    A = np.loadtxt(filename, delimiter=delim)
    return A[:,0], A[:,1]
    

def plot_data(E, G):
    """
    Plot numerical data and fit.
    """
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    En = np.logspace(np.log(E[0]), np.log(E[-1]))
    ax[0].loglog(E, G, 'b.-')
    for i in range(len(E)):
        ax[0].text(E[i], G[i], str(i), horizontalalignment='center', verticalalignment='bottom').set_clip_on(True)
    
    ax[0].set_xlabel('Photon energy (MeV)')
    ax[0].set_ylabel('Flux')

    ymx = np.ceil(np.log10(np.amax(G)))
    ymn = np.floor(np.log10(np.amin(G)))
    xmx = np.ceil(np.log10(np.amax(E)))
    xmn = np.floor(np.log10(np.amin(E)))
    ax[0].set_ylim([10**(ymn), 10**ymx])
    ax[0].set_xlim([10**xmn, 10**xmx])

    ax[1].semilogx(E, G, 'b.-')
    for i in range(len(E)):
        ax[1].text(E[i], G[i], str(i), horizontalalignment='center', verticalalignment='bottom').set_clip_on(True)

    ax[1].set_xlabel('Photon energy (MeV)')
    ax[1].set_ylabel('Flux')

    ymx = np.ceil(np.log10(np.amax(G)))
    xmx = np.ceil(np.log10(np.amax(E)))
    xmn = np.floor(np.log10(np.amin(E)))
    ax[1].set_ylim([0, 10**ymx])
    ax[1].set_xlim([10**xmn, 10**xmx])

    fig.tight_layout()
    plt.show()


def plot_spectrum(E, G, phi, c1, c2, c3, I):
    """
    Plot numerical data and fit.
    """
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    En = np.logspace(np.log(E[0]), np.log(E[-1]))
    ax[0].loglog(E, G, 'k-.')
    ax[0].loglog(En, np.exp(func(En, phi, c1, c2, c3)), 'r-')

    ax[0].set_xlabel('Photon energy (MeV)')
    ax[0].set_ylabel('Flux')

    ymx = np.ceil(np.log10(np.amax(G)))
    ymn = np.floor(np.log10(np.amin(G)))
    xmx = np.ceil(np.log10(np.amax(E)))
    xmn = np.floor(np.log10(np.amin(E)))
    ax[0].set_ylim([10**(ymn), 10**ymx])
    ax[0].set_xlim([10**xmn, 10**xmx])

    ax[1].semilogx(E, G, 'k-.')
    ax[1].semilogx(En, np.exp(func(En, phi, c1, c2, c3)), 'r-')

    ax[1].set_xlabel('Photon energy (MeV)')
    ax[1].set_ylabel('Flux')

    ymx = np.ceil(np.log10(np.amax(G)))
    xmx = np.ceil(np.log10(np.amax(E)))
    xmn = np.floor(np.log10(np.amin(E)))
    ax[1].set_ylim([0, 10**ymx])
    ax[1].set_xlim([10**xmn, 10**xmx])

    fig.suptitle(f'phi={np.exp(phi)*I:.3e}, c1={c1:.3f}, c2={c2:.3f}, c3={c3:.3f}')
    fig.tight_layout()
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Fit function to numerical Compton photon spectrum.

Fits a Compton photon spectrum to the fit used in DREAM. The spectrum should be
given as a function of photon energy, in MeV units. The data should be stored in
a CSV or TSV file, with energies in the first column and photon fluxes in the
second column. By default, the script assumes that photon fluxes are given as
'photons per energy bin' and adjusts the fitting accordingly."""
    )

    parser.add_argument('spectrum', help="File containing photon spectrum data.")
    parser.add_argument('-n', '--no-normalize-per-bin', help="Do not normalize the energy spectrum to the energy bin size.", action='store_true')
    parser.add_argument('-p', '--plot', help="Plot the fit and numerical data.", action='store_true')
    parser.add_argument('-s', '--scale', help="Factor by which to multiply the input data.", nargs='?', default=1, type=float)
    parser.add_argument('--plot-data', help="Plot the numerical data with indices.", action='store_true')
    parser.add_argument('--skip', help="Skip elements with the specified indices.", nargs='*', type=int, default=[])
    parser.add_argument('--skip-below', help="Skip all values below this index.", nargs='?', type=int, default=0)
    parser.add_argument('--skip-above', help="Skip all values above this index.", nargs='?', type=int, default=-1)

    return parser.parse_args()


def main():
    args = parse_args()

    if args.spectrum[-4:] == '.csv':
        E, G = load_csv(args.spectrum)
    elif args.spectrum[-4:] == '.tsv':
        E, G = load_tsv(args.spectrum)
    else:
        raise Exception("Unrecognized file type of spectrum file.")

    if args.plot_data:
        plot_data(E, G)

    # Normalize the spectrum to the bin size
    if not args.no_normalize_per_bin:
        dE = np.zeros(E.size)
        dE[0] = E[0]
        dE[1:] = E[1:] - E[0:-1]
        G /= dE

    # Remove data points?
    E = E[args.skip_below:]
    G = G[args.skip_below:]

    if args.skip_above >= 0:
        E = E[:args.skip_above]
        G = G[:args.skip_above]

    if args.skip:
        # Transform indices after 'skip_below' and 'skip_above'
        l = []
        for i in args.skip:
            if i >= args.skip_below:
                if args.skip_above >= 0 and i > args.skip_above:
                    continue
                else:
                    l.append(i - args.skip_below)

        print(f'Removing indices {l}')
        E = np.delete(E, l)
        G = np.delete(G, l)

    G *= args.scale

    params = fit_spectrum(E, G)
    I = eval_integral(*params) / np.exp(params[0])

    print('Best fitting parameters:')
    print(f'  PHI = {np.exp(params[0])*I:.3e}')
    print(f'   C1 = {params[1]:.3f}')
    print(f'   C2 = {params[2]:.3f}')
    print(f'   C3 = {params[3]:.3f}')

    if args.plot:
        plot_spectrum(E, G, *params, I=I)

    return 0


if __name__ == '__main__':
    sys.exit(main())


