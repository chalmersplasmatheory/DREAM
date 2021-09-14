#!/usr/bin/env python3
#
# FIT KINETIC IONIZATION CROSS SECTION TO ADAS
# --------
# This script tries to fit the total electron impact ionization cross-section
# (ICS) so that in the limit of a Maxwellian electron distribution function, it
# agrees exactly with the ADAS ionization coefficients (SCD).
#
# The fitting is done by modifying the ICS to include 2-4 free parameters. The
# parameters are then adjusted using a non-linear least-squares algorithm until
# the integral of sigma*fMe, where fMe is a Maxwellian, agrees with the
# available ADAS data at all temperatures.
#
# (This script was translated and automated by Mathias Hoppe from a Matlab
# version originally written by Ola Embreus)
# #################################################################################

import ADAS
import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import sys


PARAMETERS = {}


def nearest(arr, val):
    """
    Locate the element in the given array 'arr' which is
    closest to the specified value 'val'.
    """
    return arr[np.abs(arr-val).argmin()]


def fitSpecies(species, method='single', visualize=False):
    """
    Fits the total electron impact ICS to the ADAS coefficients for the
    given species.

    :param str species:    Name of atomic species to fit ICS for.
    :param str method:     Name of method to use for fitting: 'single', 'single_3p', or 'double'.
    :param bool visualize: If ``True``, plots the fit against ADAS data for each charge state.
    """
    I_scd, Z, _, T, _ = ADAS.data.getIonizationData(species)

    # Fit charge states Z0 = 0, 1, ..., Z-1
    fit = []
    for Z0 in range(Z):
        if Z0 == 0:
            print('{}'.format(Z0), end="", flush=True)
        else:
            print(', {}'.format(Z0), end="", flush=True)

        if species in PARAMETERS and str(Z0) in PARAMETERS[species]:
            p = PARAMETERS[species][str(Z0)]
            Tl, Tu = nearest(T, p['Tlow']), nearest(T, p['Tupp'])

            _, _, output = ADAS.fitKineticIonizationForSpecies(species, Z0, method, T_lower=Tl, T_upper=Tu)
        else:
            print("WARNING: No temperature limits set for species '{}' in charge state Z0={}. This may lead to poor fitting results.".format(species, Z0))
            print("(use the script 'ADAS/gui/ManualFit.py' to select appropriate temperature limits and add to 'kineticCrossSectionParams.json')")

            _, _, output = ADAS.fitKineticIonizationForSpecies(species, Z0, method)

        fit.append(output)

        # Show fit?
        if visualize:
            I_fit = ADAS.fit.evaluateAveragedCrossSection(method=method, T=T, **output)

            plt.loglog(T, I_scd[Z0,:,0], 'k', label='ADAS')
            plt.loglog(T, I_fit, 'r--', label='Fit')
            plt.loglog(T, I_fit, 'rx')
            plt.title('{} -- $Z_0 = {}$'.format(species, Z0))
            plt.legend()

            plt.show()

    print('.')

    return fit, Z


def compile_fits(fits, Z, method, outputfile, inttype='len_t', realtype='real_t'):
    """
    Compile fitted kinetic ionization cross-sections to C++ files.
    """
    s  = "/**\n"
    s += " * This file was automatically generated by 'tools/fit_kineticionization.py'.\n"
    s += " */\n\n"

    s += "/* Number of parameters used in fit */\n"
    if method == 'single':
        s += "const len_t nParams = 3;\n\n"
    elif method == 'single3p':
        s += "const len_t nParams = 4;\n\n"
    elif method == 'double':
        s += "const len_t nParams = 6;\n\n"
    else:
        raise Exception("Unrecognized fitting method: '{}'.".format(method))

    for e in fits.keys():
        f = fits[e]

        s += "/* {} */\n".format(e)
        s += "const {} {}_Z = {};\n".format(inttype, e, Z[e])
        s += "const {0} {1}_params[nParams*{1}_Z] = {{\n".format(realtype, e)

        if method in ['single', 'single3p']:
            s += "/*               C                     I                   betaStar */\n"
        elif method == 'double':
            s += "/*               C1                  C2                    I1                     I2                  betaStar           beta2 */\n"

        for Z0 in range(len(f)):
            s += "/* Z0 = {:2d} */   ".format(Z0)

            if method in ['single', 'single3p']:
                s += "{:18.15f}, {:21.15f}, {:17.15f}".format(f[Z0]['C1'], f[Z0]['DI1'], f[Z0]['betaStar'])
            elif method == 'double':
                s += "{:18.15f}, {:18.15f}, {:21.15f}, {:21.15f}, {:17.15f}, {:17.15f}".format(f[Z0]['C1'], f[Z0]['C2'], f[Z0]['DI1'], f[Z0]['DI2'], f[Z0]['betaStar'], f[Z0]['beta2'])

            if Z0+1 < Z[e]:
                s += ","

            s += "\n"

        s += "};\n\n"

    s += "const {} IonKineticIonizationTerm::nParamsForFit = nParams;\n".format(inttype)
    s += "const {} IonKineticIonizationTerm::kinetic_rate_n = {};\n".format(inttype, len(fits))
    s += "struct IonKineticIonizationTerm::kinetic_ionization_rate IonKineticIonizationTerm::kinetic_rate_table[{}] = {{\n".format(len(fits))

    els = list(fits.keys())
    for i in range(len(fits)):
        s += "    {{\"{0}\",{0}_Z,{0}_params}}".format(els[i])

        if i+1 < len(fits):
            s += ","

        s += "\n"

    s += "};\n\n"

    if outputfile is not None:
        # Create directory if it doesn't already exist
        pathlib.Path(outputfile).parent.mkdir(parents=True, exist_ok=True)

        # Write C++ file
        with open(outputfile, 'w') as f:
            f.write(s)

    return s


def load_parameters():
    """
    Load specific fit parameters from the standard JSON file.
    """
    paramfile = (pathlib.Path(__file__).parent / 'kineticCrossSectionParams.json').resolve()
    with open(paramfile, 'r') as f:
        params = json.load(f)

    return params


def main():
    """
    Program entry point.
    """
    path = pathlib.Path(__file__).parent.resolve()
    elementsfile = str((path / 'elements.json').resolve())
    outputfile = str((path / '../src/Atomics/kineticionizationdata.cpp').resolve())

    parser = argparse.ArgumentParser(description="Fit kinetic ionization cross-section to ADAS ionization coefficients")
    parser.add_argument('-e', '--elements', dest='elements', action='extend', nargs='*', help="Name of elements to include in fitted kinetic ionization cross-section data.")
    parser.add_argument('--elementsfile', dest='elementsfile', action='store', default=elementsfile, type=str, help="Name of file containing ADAS element specifications.")
    parser.add_argument('--method', dest='method', action='store', type=str, default='single', help="Fitting method to use (number of free parameters): 'single', 'single3p' or 'double'.")
    parser.add_argument('-o', '--output', dest='output', action='store', default=outputfile, help="Name of output C++ source file to generate.")
    parser.add_argument('--type-int', dest='inttype', action='store', default='len_t', help="C++ type to use for integers.")
    parser.add_argument('--type-real', dest='realtype', action='store', default='real_t', help="C++ type to use for real numbers.")

    args = parser.parse_args()

    if args.elements is not None:
        els = [x.lower().capitalize() for x in args.elements]
        print('Fitting ',end="")
        for i in range(len(els)):
            if i > 0:
                print(', {}'.format(els[i]), end="")
            else:
                print('{}'.format(els[i]), end="")

        print("\n")
    else:
        els = []

        e = list(ADAS.data.ELEMENTS.keys())
        print('Fitting ',end="")
        for i in range(len(ADAS.data.ELEMENTS.keys())):
            if i > 0:
                print(', {}'.format(e[i]), end="")
            else:
                print('{}'.format(e[i]), end="")

        print("\n")

    # Fit each charge state of each species...
    fits, Z = {}, {}
    for e in  ADAS.data.ELEMENTS.keys():
        # Skip elements not requested by user...
        # (but only if user requested > 0 elements)
        if len(els) > 0:
            if e.lower().capitalize() not in els:
                continue

        print('{:2s}: '.format(e), end="", flush=True)
        fits[e], Z[e] = fitSpecies(e, method=args.method)

    compile_fits(fits, Z=Z, method=args.method, outputfile=args.output, inttype=args.inttype, realtype=args.realtype)

    return 0


if __name__ == '__main__':
    PARAMETERS = load_parameters()
    sys.exit(main())
    #fitSpecies('Ar', visualize=True)


