# Routines for working with ADAS data

import datetime
import json
import numpy as np
import os
import pathlib
import sys
import time
import urllib.request

sys.path.append('..')

import get_nist


CACHEDIR = None
ELEMENTS = None


def getIonizationData(species):
    """
    Get ionization data for the named ion species.
    """
    global CACHEDIR

    Z, n, T, I_scd = _loadADAS(species, 'scd')
    _, _, E_ion    = get_nist.load_elements([species], 'ionization', cache=True, cachedir=CACHEDIR)

    n     = np.power(10, n)
    T     = np.power(10, T)
    I_scd = np.power(10, I_scd)
    E_ion = E_ion[0]

    return I_scd, Z, n, T, E_ion


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


def _initADAS():
    """
    Initialize this module.
    """
    global CACHEDIR, ELEMENTS

    p = pathlib.Path(__file__).parent.parent.absolute()

    CACHEDIR = '{}/cache'.format(p)
    ELEMENTS = load_element_list(str(p / 'elements.json'))


def _loadADAS(species, datatype, cache=True):
    """
    Load the specified data for the given species from ADAS.

    :param str species:  Name of ion species to load data for.
    :param str datatype: Name of datatype to load ('acd', 'scd', 'plt', 'prb', 'ccd')

    :returns: Atomic number ``Z``, density vector ``n``, temperature vector ``T``, coefficient values ``v``
    Note that ``n``, ``T`` and ``v`` contain base-10 logarithm values and must be exponentiated before use.
    """
    global CACHEDIR, ELEMENTS

    year = ELEMENTS[species]
    if type(year) == dict:
        y = year[datatype]
    else:
        y = year

    Z, n, T, v = parse_adas(download_adas(species, y, datatype, cache=cache, cachedir=CACHEDIR))

    return Z, n, T, v


def download_adas(element, year, datatype, cache=False, cachedir=None):
    """
    Downloads data of the specified type for the specified element,
    for the given year. The 'datatype' parameter may be either of the following:

      acd  -- Effective recombination coefficients
      scd  -- Effective ionization coefficients
      ccd  -- Charge exchange effective recombination coefficients
      plt  -- Line power driven by excitation of dominant ions
      prb  -- Continuum and line power driven by recombination
              and Bremsstrahlung of dominant ions

    If 'cache' is True, the downloaded data is stored in a text
    file in the specified cache directory. Alternatively, if the
    file already exists in the specified cache directory, data is
    read from it.
    """
    # The year may also be the name of a specific file to load
    yearIsFile = (len(year) > 3 and year[-4:] == '.dat')
    if yearIsFile:
        fname = year
    else:
        fname = '{0}{1}_{2}.dat'.format(datatype.lower(), year, element.lower())

    # Construct ADAS data url
    dt = datatype.lower()
    url = 'https://open.adas.ac.uk/download/adf11/{0}{1}/{2}'.format(dt, year, fname)

    data = None
    fpath = pathlib.PurePath(cachedir, fname)
    if (cache and os.path.isfile(fpath)) or yearIsFile:
        with open(fpath, 'r') as f:
            data = f.read()
    else:   # Load from open.adas.ac.uk
        with urllib.request.urlopen(url) as f:
            if f.status != 200:
                raise Exception("Failed to download '{}' from Open-ADAS.".format(url))

            data = f.read().decode('ascii')

        # Save data to disk?
        if cache:
            # Create cache directory if it doesn't exist
            pathlib.Path(cachedir).mkdir(parents=True, exist_ok=True)

            with open(fpath, 'w') as f:
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
        # +6: convert from cm^-3 to m^-3
        n += [float(x)+6 for x in lines[i].split()]
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
            # -6: convert from cm^3 to m^3
            ldata += [float(x)-6 for x in lines[i].split()]
            i += 1

        data.append(ldata)
        # Skip ' --'...
        i += 1

    data = np.reshape(np.array(data), (Z, nT, nn))
    n = np.array(n)
    T = np.array(T)

    return Z, n, T, data


def load_element(element, year, cache=True, cachedir=None):
    """
    Load all rate coefficient data for the specified element.
    """
    data = {'acd': None, 'ccd': None, 'scd': None, 'plt': None, 'prb': None}

    print("Loading data for element '{}'... ".format(element), end="")
    Z = 0
    t = time.time()
    for key in data:
        if type(year) == dict:
            y = year[key]
        else:
            y = year
            
        # For rates other than CCD (charge-exchange), we use
        # data for H also for its isotopes D and T.
        el = element
        if key != 'ccd' and el in ['D', 'T']:
            el = 'H'

        Z, n, T, v = parse_adas(download_adas(el, y, key, cache=cache, cachedir=cachedir))

        data[key] = {
            'n': n,
            'T': T,
            'data': v
        }

    data['Z'] = Z

    # Handle hydrogen isotope mass numbers specifically...
    if element == 'H': data['A'] = 1
    elif element == 'D': data['A'] = 2
    elif element == 'T': data['A'] = 3
    else: data['A'] = 2*Z

    print("{} ms".format((time.time()-t)*1e3))

    return data


# Initialize this module
_initADAS()

