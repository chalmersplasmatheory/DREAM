#!/usr/bin/python3 -i

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import sys

# TODO Set this up properly
path = str((pathlib.Path(__file__).parent / '..').resolve().absolute())
sys.path.append(path)

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput


def create_argparser():
    parser = argparse.ArgumentParser(description="DREAM Output CLI")

    parser.add_argument('-l', '--lazy', help="Load the output lazily and only read data on-demand", dest="lazy", action="store_true")
    parser.add_argument('-r', '--read', help="Load the output to memory immediately (no lazy read)", dest="lazy", action="store_false")
    parser.add_argument('-s', '--no-settings', help="Do not load settings from the output file", dest="settings", action="store_false")

    parser.add_argument('output', help="DREAM output file to load", type=str, nargs='?')

    parser.set_defaults(lazy=True, output='output.h5', settings=True)

    return parser.parse_args()
    

def main():
    args = create_argparser()

    do = DREAMOutput(args.output, lazy=args.lazy, loadsettings=args.settings)
    """
    if len(sys.argv) == 1:
        do = DREAMOutput('output.h5')
    elif len(sys.argv) == 2:
        do = DREAMOutput(sys.argv[1])
    else:
        print('ERROR: Invalid command line arguments. Expected at most one argument (name of file).')
        sys.exit(1)
    """

    setup_interactive(do, glob=globals())


if __name__ == '__main__':
    main()
    
    # When main returns we continue to an
    # interactive Python session...

