#!/usr/bin/env python3
#
# DOWNLOAD DATA FROM OPEN-ADAS
#
# This script downloads atomic ionization/recombination and radiation data from
# the Open-ADAS database (https://open.adas.ac.uk/).
#

import numpy as np
import os
import sys
import time
import urllib.request


# Define which datasets to use. The value indicates which year
# the dataset corresponds to. Please check the Open_ADAS documentation
# (https://open.adas.ac.uk/man/appxa-11.pdf) for the quality of
# the dataset before adding it to this list.
ELEMENTS = {
    'H': '96',
    'He': '96',
    'Li': '96',
    'Be': '96',
    'C':  '96',
    'Ne': '96',
    'Ar': '89'
}


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
            data = f.read().decode('ascii')

        # Save data to disk?
        if cache:
            with open(fname, 'w') as f:
                f.write(data)

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
    while not lines[i].startswith(' --'):
        T += [float(x) for x in lines[i].split()]
        i += 1

    # Skip ' ---' line...
    i += 1

    # Load data for each charge state
    data = []
    while lines[i][0] != 'C':
        ldata = []
        # Load data corresponding to single charge state
        while not lines[i].startswith(' --') and lines[i] != ' ':
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

    return n, T, data


def main(argv):
    data = download_adas('Ar', 89, 'SCD', cache=True)

    t = time.time()
    parse_adas(data)
    print('Duration: {} ms'.format((time.time()-t)*1e3))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

