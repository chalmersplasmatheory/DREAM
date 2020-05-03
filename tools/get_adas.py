#!/usr/bin/env python3
#
# DOWNLOAD DATA FROM OPEN-ADAS
#
# This script downloads atomic ionization/recombination and radiation data from
# the Open-ADAS database (https://open.adas.ac.uk/).
#

import argparse
import datetime
import json
import numpy as np
import os
import pathlib
import sys
import time
import urllib.request


def download_adas(element, year, datatype, cache=False):
    """
    Downloads data of the specified type for the specified element,
    for the given year. The 'datatype' parameter may be either of the following:

      acd  -- Effective recombination coefficients
      scd  -- Effective ionization coefficients
      plt  -- Line power driven by excitation of dominant ions
      prb  -- Continuum and line power driven by recombination
              and Bremsstrahlung of dominant ions

    If 'cache' is True, the downloaded data is stored in a text
    file in the current working directory. Alternatively, if the
    file already exists in the current working directory, data is
    read from it.
    """
    # Construct ADAS data url
    dt = datatype.lower()
    fname = '{0}{1}_{2}.dat'.format(datatype.lower(), year, element.lower())
    url = 'https://open.adas.ac.uk/download/adf11/{0}{1}/{2}'.format(dt, year, fname)

    data = None
    if cache and os.path.isfile(fname):
        with open(fname, 'r') as f:
            data = f.read()
    else:   # Load from open.adas.ac.uk
        with urllib.request.urlopen(url) as f:
            if f.status != 200:
                raise Exception("Failed to download '{}' from Open-ADAS.".format(url))

            data = f.read().decode('ascii')

        # Save data to disk?
        if cache:
            with open(fname, 'w') as f:
                f.write(data)

    if data.startswith('<!DOCTYPE'):
        raise Exception("Failed to download '{}' from Open-ADAS. The file does not appear to exist.".format(url))

    return data


def parse_adas(data):
    """
    Parses the given ADAS data file. This function expects the
    data to be in exactly the same format as obtained through
    the Open ADAS database.
    """
    lines = data.splitlines()

    # Get dimensions
    dims = lines[0].split()

    Z  = int(dims[0])
    nn = int(dims[1])
    nT = int(dims[2])

    # Density
    i = 2
    n = []
    while len(n) < nn:
        n += [float(x) for x in lines[i].split()]
        i += 1

    # Temperature
    T = []
    while not lines[i].strip(' C').startswith('--'):
        T += [float(x) for x in lines[i].split()]
        i += 1

    # Skip ' ---' line...
    i += 1

    # Load data for each charge state
    data = []
    while lines[i][0] != 'C':
        ldata = []
        # Load data corresponding to single charge state
        while not lines[i].strip(' C').startswith('--') and lines[i] != ' ':
            if i == 3486:
                print("'{}' is {}".format(lines[i], lines[i] == ' '))

            ldata += [float(x) for x in lines[i].split()]
            i += 1

        data.append(ldata)
        # Skip ' --'...
        i += 1

    data = np.reshape(np.array(data), (Z, nn, nT))
    n = np.array(n)
    T = np.array(T)

    return Z, n, T, data


def load_element(element, year, cache=True):
    """
    Load all rate coefficient data for the specified element.
    """
    data = {'acd': None, 'scd': None, 'plt': None, 'prb': None}

    print("Loading data for element '{}'... ".format(element), end="")
    Z = 0
    t = time.time()
    for key in data:
        Z, n, T, v = parse_adas(download_adas(element, year, key, cache=cache))

        data[key] = {
            'n': n,
            'T': T,
            'data': v
        }

    data['Z'] = Z

    print("{} ms".format((time.time()-t)*1e3))

    return data


def load_element_list(filename):
    """
    Load the list of elements in the ADAS database as
    defined for this script in the 'elemenets.json' file,
    located in the same directory as this script.

    elements.json:
      Define which datasets to use. The value indicates which year
      the dataset corresponds to. Please check the Open_ADAS documentation
      (https://open.adas.ac.uk/man/appxa-11.pdf) for the quality of
      the dataset before adding it to this list.
    """

    ELEMENTS = None
    with open(filename, 'r') as f:
        ELEMENTS = json.load(f)

    return ELEMENTS


def compile_elements(elements, outputfile=None):
    """
    Generate a C++ file containing the given ADAS data.
    """
    ds = ""
    ss = "const len_t adas_rate_n = {0};\nstruct adas_rate adas_rate_table[{0}] = {{\n".format(len(elements))
    sd = None

    # Write data
    for elname, ratedata in elements.items():
        Z = ratedata['Z']

        if sd is not None: sd += "},\n"
        else: sd = ""

        sd += "\t{{\"{0}\",{0}_Z".format(elname)

        ds += "/* {} */\n".format(elname)
        ds += "const len_t {0}_Z = {1};\n".format(elname, Z)

        # acd, scd, plt, prb
        for dt, data in ratedata.items():
            if dt == 'Z': continue

            nn = len(data['n'])
            nT = len(data['T'])

            ds += "const len_t {0}_{1}_nn = {2};\nconst len_t {0}_{1}_nT = {3};\n".format(elname, dt, nn, nT)

            ds += "const real_t {0}_{1}_n[{2}] = {{".format(elname, dt, nn)
            ds += ','.join(['{:.5f}'.format(x) for x in data['n']])
            ds += "};\n"

            ds += "const real_t {0}_{1}_T[{2}] = {{".format(elname, dt, nT)
            ds += ','.join(['{:.5f}'.format(x) for x in data['T']])
            ds += "};\n"

            ds += "const real_t {0}_{1}_coeff[{2}] = {{".format(elname, dt, Z*nn*nT)
            ds += ','.join(['{:.5f}'.format(x) for x in data['data'].flatten()])
            ds += "};\n\n"

            sd += ",{0}_{1}_nn,{0}_{1}_nT,{0}_{1}_n,{0}_{1}_T,{0}_{1}_coeff".format(elname, dt)

    sd += "}\n};\n"

    filecontents  = "/* This file was auto-generated by 'get_adas.py' on {} */\n\n".format(datetime.datetime.now().isoformat(sep=' ', timespec='seconds'))
    filecontents += "#include \"DREAM/adasdata.h\"\n\n"
    filecontents += ds + ss + sd

    if outputfile is not None:
        # Create directory if it doesn't already exist
        pathlib.Path(outputfile).parent.mkdir(parents=True, exist_ok=True)

        # Write C++ file
        with open(outputfile, 'w') as f:
            f.write(filecontents)

    return filecontents


def main(argv):
    path = pathlib.Path(__file__).parent.absolute()
    elementsfile = '{}/elements.json'.format(path)
    outputfile = os.path.abspath('{}/../src/ADAS/adasdata.cpp'.format(path))

    parser = argparse.ArgumentParser(description="Download and compile rate coefficients from Open-ADAS.")
    parser.add_argument('--elements', dest='elementsfile', action='store', default=elementsfile, type=str, help="Name of file containing ADAS element specifications.")
    parser.add_argument('--no-cache', dest='cache', action='store_false', help="Forces data to be downloaded from Open-ADAS and prevents Open-ADAS files from being stored locally.")
    parser.add_argument('--no-compile', dest='compile', action='store_false', help="Do not generate C++ source files with the rate coefficients.")
    parser.add_argument('-o', '--output', dest='output', action='store', default=outputfile, help="Name of output C++ source file to generate.")

    args = parser.parse_args()

    ADAS_ELEMENT_LIST = load_element_list(args.elementsfile)

    ELEMENTS = {}
    for element, year in ADAS_ELEMENT_LIST.items():
        ELEMENTS[element] = load_element(element, year, cache=args.cache)

    # Compile 
    if not args.compile:
        return 0

    compile_elements(ELEMENTS, outputfile=args.output)
        
    #data = download_adas('Ar', 85, 'SCD', cache=True)

    #t = time.time()
    #parse_adas(data)
    #print('Duration: {} ms'.format((time.time()-t)*1e3))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

