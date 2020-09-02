#!/usr/bin/env python3

import argparse
import pathlib
import sys
import time

import dreamtests

try:
    import dreamtests
except:
    sys.path.append(pathlib.Path(__file__).parent.absolute())
    import dreamtests

try:
    import DREAM
except ImportError:
    sys.path.append(str((pathlib.Path(__file__).parent / '..' / '..' / 'py').resolve().absolute()))
    import DREAM


# Import test modules
from code_conductivity import code_conductivity
from code_runaway import code_runaway
from ts_adaptive import ts_adaptive
from DREAM_avalanche import DREAM_avalanche


TESTS = [
    'code_conductivity',
    'code_runaway',
    'DREAM_avalanche',
    'ts_adaptive'
]


def print_help():
    """
    Prints command help.
    """
    global TESTS

    print("Physics tests for DREAM\n")

    print("Usage:")
    print("    runtests.py           Show this help message.")
    print("    runtests.py all       Run all tests.")
    print("    runtests.py [FLAGS] [test1 [test2 [...]]]")
    print("                          Run the tests with names 'test1', 'test2' etc.\n")

    print("Options:")
    print("    --plot                Plot result instead of comparing automatically\n")

    print("Available tests:")

    for test in TESTS:
        print("    {}".format(test))


def runtest(name, args):
    """
    Run the named test.
    """
    print("\x1B[1m:: {} \x1B[0m".format(name))

    if name not in globals():
        print("ERROR: Unrecognized test: '{}'".format(name))
        return

    globals()[name].run(args)


def runall(args):
    """
    Run all available tests.
    """
    global TESTS

    for test in TESTS:
        runtest(test, args)

    
def main(argv):
    """
    Program entry point.
    """
    parser = argparse.ArgumentParser(description='DREAM physics tests')

    parser.add_argument('--plot', help="In tests where applicable, plot results instead of comparing automatically", action="store_true")
    parser.add_argument('tests', help="List of tests to run", type=str, nargs='*')

    args = parser.parse_args()
    arglist = {'plot': args.plot}

    if len(args.tests) == 0:
        print_help()
        return 1
    elif len(args.tests) == 1 and args.tests[0].lower() == 'all':
        runall(arglist)
    else:
        for test in args.tests:
            runtest(test, arglist)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


