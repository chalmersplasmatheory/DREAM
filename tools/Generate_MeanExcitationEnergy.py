#!/usr/bin/env python3
#
# GENERATE MEAN EXCITATION ENERGY DATA FOR Z>18
# --------
# This script generates Mean Excitation Energy (MEE) data from
# semi-empirical approximation for elements with Z>18 and store them in C++ file.
# This is used in DREAM to complement data from Sauer et al.
#
# J. Walkowiak
# 08.2022
# #################################################################################

import numpy as np
import scipy.integrate as integrate
import math
import sys
import pathlib
import argparse

def GenerateMEE(Z_min=19, Z_max=86):
    P = np.array([[1.1831,	0.1738,	0.0913,	0.0182,	0.7702],
        [0.8368,	1.0987,	0.9642,	1.2535,	0.2618],
        [0.3841,	0.6170,	1.0000,	1.0000,	1.0000],
        [0.5883,	0.0461,	1.0000,	1.0000,	1.0000]])

    lnI = np.empty((Z_max+1, Z_max+1)) * np.nan
    R = np.logspace(-4, 1, 10000)

    for Z in range(Z_min, Z_max+1):
        lambda_Z =  np.array([ P[0,0]*Z**P[1,0], P[0,1]*Z**P[1,1], P[0,2]*Z**P[1,2], P[0,3]*Z**P[1,3], P[0,4]*Z**P[1,4]])
        n_Z =       np.array([ P[2,0]*Z**P[3,0], P[2,1]*Z**P[3,1], P[2,2]*Z**P[3,2], P[2,3]*Z**P[3,3], P[2,4]*Z**P[3,4]])
        for N in range(1, Z+1):
            q = (Z-N) / Z
            lambda_i = lambda_Z * ( (1 - q**(n_Z+1)) / (1-q) ) ** (1/2)

            e_density = (lambda_i[0] ** 2 * min(N, 2) * np.exp(-lambda_i[0] * R))   \
                        + (lambda_i[1] ** 2 * min( max(N-2,0) ,8) * np.exp(-lambda_i[1] * R))   \
                        + (lambda_i[2] ** 2 * min( max(N-10,0) ,18) * np.exp(-lambda_i[2] * R)) \
                        + (lambda_i[3] ** 2 * min( max(N-28,0) ,26) * np.exp(-lambda_i[3] * R)) \
                        + (lambda_i[4] ** 2 * max(N-54,0) * np.exp(-lambda_i[4] * R))

            e_density = 1 / (4*math.pi*R) * e_density

            local_frequency = np.log(np.sqrt(4*math.pi*e_density), out=np.zeros_like(R), where=e_density!=0)

            lnI_au = 4 * math.pi / N * integrate.trapezoid(R**2 * e_density * local_frequency, R) \
                     + math.exp(1-(N/Z)) - 0.9
            lnI[Z,N] = lnI_au + math.log(27.211)

    MEE = np.exp(lnI)

    return MEE


def element_symbol(Z):
    """
    Return element symbol by atomic number Z
    """
    elements_list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                     "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
                     "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                     "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
                     "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
                     "Po", "At", "Ra"]
    return elements_list[Z-1]


def saveMEE2cpp(MEE, outputfile, Z_max=86, realtype='real_t'):
    """
    Save data to C++ file.
    """
    s  = "/**\n"
    s += " * This file was automatically generated by 'tools/Generate_MeanExcitationEnergy.py'.\n"
    s += " */\n\n"

    s += "/*  List of mean excitation energies in units of eV  */\n"
    s += "const {0} MEAN_EXCITATION_ENERGY_EXTENDED[{1}][{1}] = {{\n".format(realtype, Z_max)

    for Z in range(1, Z_max+1):

        s += "/* {} */ {{".format(element_symbol(Z))
        for Z0 in range(1, Z_max+1):
            s += " {:.1F}".format(MEE[Z][Z0])
            if Z0<Z_max: s += ','

        if Z < Z_max:   s += ' },\n'
        else:           s += ' }\n'
    s += "};\n\n"

    if outputfile is not None:
        # Create directory if it doesn't already exist
        pathlib.Path(outputfile).parent.mkdir(parents=True, exist_ok=True)

        # Write C++ file
        with open(outputfile, 'w') as f:
            f.write(s)

    return s


def main():
    """
    Program entry point.
    """
    path = pathlib.Path(__file__).parent.resolve()
    outputfile = str((path / '../src/Equations/MeanExcitationEnergy_Extended.cpp').resolve())
    #outputfile = str((path / './MeanExcitationEnergy_Extended.cpp').resolve())

    parser = argparse.ArgumentParser(description="Generate Mean Excitation Energy (MEE) data from semi-empirical approximation")
    parser.add_argument('-Zmin', '--Z_min', dest='Z_min', action='store', default=19, help="Minimum atomic number for which MEE will be generated")
    parser.add_argument('-Zmax', '--Z_max', dest='Z_max', action='store', default=86, help="Maximum atomic number for which MEE will be generated (max 86)")
    parser.add_argument('-o', '--output', dest='output', action='store', default=outputfile, help="Name of output C++ source file to generate.")
    parser.add_argument('--type-real', dest='realtype', action='store', default='real_t', help="C++ type to use for real numbers.")

    args = parser.parse_args()

    saveMEE2cpp(GenerateMEE(Z_min=args.Z_min, Z_max=args.Z_max), outputfile=args.output, Z_max=args.Z_max, realtype=args.realtype)

    return 0


if __name__ == '__main__':
    sys.exit(main())


