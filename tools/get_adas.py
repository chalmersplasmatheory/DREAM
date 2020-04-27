#!/usr/bin/env python3
#
# DOWNLOAD DATA FROM OPEN-ADAS
#
# This script downloads atomic ionization/recombination and radiation data from
# the Open-ADAS database (https://open.adas.ac.uk/).
#

import numpy as np
import sys
import time
import urllib.request


# Define which datasets to use
ELEMENTS_SCD = {
    'H': '93',
    'He': '96',
    'Li': '96',
    'Be': '93',
    'Ar': '89'
}


def download_adas(element, year, datatype):
    """
    Downloads data of the specified type for the specified element,
    for the given year. The 'datatype' parameter may be either of the following:

      acd  -- Effective recombination coefficients
      scd  -- Effective ionization coefficients
      plt  -- Line power driven by excitation of dominant ions
      prb  -- Continuum and line power driven by recombination
              and Bremsstrahlung of dominant ions
    """
    # Construct ADAS data url
    dt = datatype.lower()
    fname = '{0}{1}_{2}.dat'.format(datatype.lower(), year, element.lower())
    url = 'https://open.adas.ac.uk/download/adf11/{0}{1}/{2}'.format(dt, year, fname)

    data = None
    with urllib.request.urlopen(url) as f:
        data = f.read().decode('ascii')

    # DEBUG Temporarily store data on disk
    with open(fname, 'w') as f:
        f.write(data)


def parse_adas(data):
    """
    Parses the given ADAS data.
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
    #download_adas('Ar', 89, 'SCD')
    # DEBUG Load in pre-saved data
    data = None
    with open('scd89_ar.dat', 'r') as f:
        data = f.read()

    t = time.time()
    parse_adas(data)
    print('Duration: {} ms'.format((time.time()-t)*1e3))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

