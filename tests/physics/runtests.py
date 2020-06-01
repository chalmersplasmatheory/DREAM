#!/usr/bin/env python3

import pathlib
import sys
import time

"""
try:
    import dreamtests
except:
    sys.path.append(pathlib.Path(__file__).parent.absolute())
    import dreamtests
"""

try:
    import DREAM
except ImportError:
    sys.path.append(str((pathlib.Path(__file__).parent / '..' / '..' / 'py').resolve().absolute()))
    import DREAM


# Import test modules
from code_conductivity import code_conductivity


TESTS = [
    'code_conductivity'
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
    print("    runtests.py [test1 [test2 [...]]]")
    print("                          Run the tests with names 'test1', 'test2' etc.\n")

    print("Available tests:")

    for test in TESTS:
        print("    {}".format(test))


def print_ok(msg):
    """
    Prints the given message with an [OK] prefix.
    """
    print("\x1B[1;32m[OK]\x1B[0m    --> {}".format(msg))


def runtest(name):
    """
    Run the named test.
    """
    print("\x1B[1m:: {} \x1B[0m".format(name))

    globals()[name].run()


def runall():
    """
    Run all available tests.
    """
    global TESTS

    for test in TESTS:
        runtest(test)

    
def main(argv):
    """
    Program entry point.
    """
    if len(argv) == 0:
        print_help()
        return 1
    elif len(argv) == 1 and argv[0].lower() == 'all':
        runall()
    else:
        for test in argv:
            runtest(test)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


