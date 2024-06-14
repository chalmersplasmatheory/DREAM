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
    popt, _ = curve_fit(func, E, np.log(G), p0=(np.amax(G), 1.2, 0.8, 0)) 
    return popt


def load_csv(filename): return load_delim(filename, ',')
def load_tsv(filename): return load_delim(filename, '\t')


def load_delim(filename, delim):
    """
    Load photon spectrum from a delimited text file.
    """
    A = np.loadtxt(filename, delimiter=delim)
    return A[:,0], A[:,1]
    

def plot_spectrum(E, G, phi, c1, c2, c3, I):
    """
    Plot numerical data and fit.
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    En = np.logspace(np.log(E[0]), np.log(E[-1]))
    ax.loglog(E, G, 'k-.')
    ax.loglog(En, np.exp(func(En, phi, c1, c2, c3)), 'r-')

    ax.set_xlabel('Photon energy (MeV)')
    ax.set_ylabel('Flux')

    mx = np.ceil(np.log10(np.amax(G)))
    ax.set_ylim([10**(mx-10), 10**mx])

    ax.set_title(f'phi={np.exp(phi)*I:.3e}, c1={c1:.3f}, c2={c2:.3f}, c3={c3:.3f}')

    fig.tight_layout()
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser('Fit function to numerical Compton photon spectrum.')

    parser.add_argument('spectrum', help="File containing photon spectrum data.")
    parser.add_argument('-p', '--plot', help="Plot the fit and numerical data.", action='store_true')
    parser.add_argument('-s', '--scale', help="Factor by which to multiply the input data.", nargs='?', default=1, type=float)

    return parser.parse_args()


def main():
    args = parse_args()

    if args.spectrum[-4:] == '.csv':
        E, G = load_csv(args.spectrum)
    elif args.spectrum[-4:] == '.tsv':
        E, G = load_tsv(args.spectrum)
    else:
        raise Exception("Unrecognized file type of spectrum file.")

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


