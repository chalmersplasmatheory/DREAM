#!/usr/bin/env python3
#
# GET NIST DATA
# This script downloads atomic data from the American National Institute of
# Standards and Technology's (NIST) Atomic Spectra Database (ASD). It generates
# a request in the manner of the Ionization Energy form at
#
#   https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
#
# and parses the response HTML.
#
# ###################

import argparse
import datetime
import json
import numpy as np
import os
import pathlib
import sys
import time
import urllib.request
import ADAS.io


def download_nist(elements, datatype='binding', cache=False, cachedir=None):
    """
    Downloads the HTML result of making a request using the ionization energy
    data form from the NIST ASD.

    elements: List of elements to download data for. This should be the name of
              the atom to download data for (e.g. H, He, Ne, Ar etc.)
    datatype: Type of data to download (either 'binding' for total binding
              energy, or 'ionization' for ionization energy)
    """
    dt = '0' if datatype == 'ionization' else '1'
    data = {
        'encodelist': 'XXT2',           # Hard-coded value
        'spectra': ','.join(elements),  # List of elements to retrieve
        'units': '1',                   # Energy units (1 = eV, 0 = cm-1, 2 = Rydberg
        'format': '1',                  # Output format (0 = HTML, 1 = ASCII)
        'order': '0',                   # Ordering (0 = by Z, 1 = by sequence)
        'at_num_out': 'on',             # We want the atomic charge number
        'ion_charge_out': 'on',         # We also want the ion charge as X data
        'sp_name_out': 'on',            # Spectrum name
        'e_out': dt,                    # 0 = ionization energy, 1 = total binding energy
        'submit': 'Retrieve+Data'       # Name of button
    }

    fname = 'nist_{}.html'.format(datatype)
    url   = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl'

    if cachedir is not None:
        fpath = str(pathlib.PurePath(cachedir, fname))
    else:
        fpath = ''

    if cache and os.path.isfile(fpath):
        with open(fpath, 'r') as f:
            nistdata = f.read()
    else:
        dd  = urllib.parse.urlencode(data).encode()
        #req = urllib.request.Request(url, data=
        with urllib.request.urlopen(url, data=dd) as f:
            if f.status != 200:
                raise Exception("Failed to download '{}' from NIST ADS.".format(url))

            nistdata = f.read().decode('ascii')

        # Save data to disk?
        if cache:
            pathlib.Path(cachedir).mkdir(parents=True, exist_ok=True)

            with open(fpath, 'w') as f:
                f.write(nistdata)

    return nistdata


def remove_html(s):
    """
    Removes all HTML tags from a string.
    """
    tmp = s
    #while (sub := tmp.find('<')) >= 0:
    sub = tmp.find('<')
    while sub >= 0:
        tmp1 = tmp[:sub]
        tmp2 = tmp[sub:]

        sub = tmp2.find('>')
        tmp = tmp1 + tmp2[sub+1:]

        sub = tmp.find('<')

    return tmp


def parse_data(nistdata):
    """
    Parses data downloaded from the NIST Atomic Spectra Database.
    The data is expected to be of type 'ASCII' and embedded in an
    HTML document.
    """
    # Retrieve '<pre></pre>' section
    preidx = nistdata.find('<pre>')
    substr = nistdata[preidx+6:]
    preidx = substr.find('</pre>')
    substr = substr[:preidx-1]

    names  = []
    Z      = []
    data   = []

    # Number of lines to remove from top of table (table header)
    N_REMOVE_TOP = 3
    # Number of lines to remove from bottom of table
    N_REMOVE_BOTTOM = 1

    lines  = substr.split("\n")
    while lines[-N_REMOVE_BOTTOM][:4] != '----':
        N_REMOVE_BOTTOM += 1

    lines  = lines[N_REMOVE_TOP:-N_REMOVE_BOTTOM]
    arr    = []
    for line in lines:
        # Separate X and Y values
        s = [x.strip() for x in line.split('|')]

        # Get charge number, state and name
        _Z    = int(s[0])
        _Z0   = int(s[2])
        _name = s[1].split()[0]

        # Filter out HTML (and surrounding ()/[])
        _v = remove_html(s[3])
        if _v[0] == '[' or _v[0] == '(':
            _v = _v[1:-1]
        _v = float(_v)

        #data.append((name, Z, Z0, v))

        if len(Z) == 0 or Z[-1] != _Z:
            # Push old data
            if arr:
                data.append(arr)

            # Add new elements
            Z.append(_Z)
            names.append(_name)
            arr = []

        arr.append(_v)

    data.append(arr)

    return (names, Z, data)


def load_elements(elements, datatype='binding', cache=False, cachedir=None):
    """
    Load data for the named elements.
    """
    data = download_nist(elements=elements, datatype=datatype, cache=cache, cachedir=cachedir)

    names, Z, data = parse_data(data)

    # Check that 'elements' is a subset of 'names'
    for e in elements:
        if e not in names:
            if cache:
                # Force download of elements
                return load_elements(elements=elements, datatype=datatype, cache=False)
            else:
                raise Exception("Unable to load NIST data for element '{}'. Element data unavailable.".format(e))

    # Check whether to return all data or a subset
    if len(elements) == len(names):
        return names, Z, data
    else:
        # Select only requested elements
        nnames, nZ, ndata = [], [], []
        for n,z,d in zip(names,Z,data):
            if n in elements:
                nnames.append(n)
                nZ.append(z)
                ndata.append(d)

        return nnames, nZ, ndata


def compile_data(nistdata, outputfile, datatype='binding', inttype='int', realtype='double'):
    """
    Generate a C++ file containing the given NIST data.
    """
    names, Z, data = nistdata
    ds = ""
    ss = "const {0} nist_{2}_n = {1};\nstruct nist_data nist_{2}_table[{1}] = {{\n".format(inttype, len(data), datatype)
    sd = None

    # Write data
    for i in range(len(data)):
        if sd is not None: sd += "},\n"
        else: sd = ""
        
        sd += "\t{{\"{0}\",{0}_{1}_Z".format(names[i], datatype)

        ds += "/* {} */\n".format(names[i])
        ds += "const {0} {1}_{2}_Z = {3};\n".format(inttype, names[i], datatype, Z[i])

        ds += "const {0} {1}_{2}_data[{3}] = {{".format(realtype, names[i], datatype, len(data[i]))
        ds += ','.join(['{:.10f}'.format(x) for x in data[i]])
        ds += "};\n"

        sd += ",{0}_{1}_data".format(names[i], datatype)

    sd += "}}\n}};\n".format(names[i], datatype)

    filecontents  = "/* This file was auto-generated by 'get_nist.py' on {} */\n\n".format(datetime.datetime.now().isoformat(sep=' ', timespec='seconds'))
    filecontents += "#include \"DREAM/nistdata.h\"\n\n"
    filecontents += ds + ss + sd

    if outputfile is not None:
        # Create directory if it doesn't already exists
        pathlib.Path(outputfile).parent.mkdir(parents=True, exist_ok=True)

        # Write C++ file
        with open(outputfile, 'w') as f:
            f.write(filecontents)

    return filecontents


def main():
    """
    Program entry point.
    """
    path       = pathlib.Path(__file__).parent.absolute()
    cachedir   = '{}/cache'.format(path)
    #elements   = ['H', 'He', 'Be', 'Ne', 'Ar']

    with open('{}/elements.json'.format(path), 'r') as f:
        elements = list(json.load(f).keys())

    if 'D' in elements: elements.remove('D')
    if 'T' in elements: elements.remove('T')

    parser = argparse.ArgumentParser(description="Download and compile ionization energies from NIST ADS")
    parser.add_argument('--cachedir', dest='cachedir', action='store', default=cachedir, type=str, help="Path to directory in which to store/load cached data files to/from.")
    parser.add_argument('--hdf5', dest='hdf5', action='store', type=str, help="Store data in the named HDF5 file")
    parser.add_argument('--ionization', dest='bindingenergy', action='store_false', help="Downloads ionization energy data instead of binding energy data")
    parser.add_argument('--no-cache', dest='cache', action='store_false', help="Forces data to be downloaded from the NIST ADS and prevents files from being stored locally.")
    parser.add_argument('--no-compile', dest='compile', action='store_false', help="Do not generate C++ source files with the ionization data.")
    parser.add_argument('-o', '--output', dest='output', action='store', default='', help="Name of output C++ source file to generate.")
    parser.add_argument('--type-int', dest='inttype', action='store', default='len_t', help="C++ type to use for integers.")
    parser.add_argument('--type-real', dest='realtype', action='store', default='real_t', help="C++ type to use for real numbers.")

    args = parser.parse_args()

    datatype = 'binding' if args.bindingenergy else 'ionization'

    data = load_elements(elements, datatype=datatype, cache=args.cache, cachedir=args.cachedir)

    if args.compile:
        if args.output == '':
            args.output = os.path.abspath('{}/../src/Atomics/nistdata_{}.cpp'.format(path, datatype))

        # Compile data to C++
        compile_data(data, outputfile=args.output, datatype=datatype, inttype=args.inttype, realtype=args.realtype)

    # Store in HDF5
    if args.hdf5:
        d= {}
        for i in range(len(data[0])):
            ion = data[0][i]
            d[ion] = {
                'Z': data[1][i],
                'data': data[2][i]
            }
        ADAS.io.save_dict(d, outputfile=args.hdf5)

    return 0


if __name__ == '__main__':
    sys.exit(main())


