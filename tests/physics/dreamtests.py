# Helper routines for the test framework

import sys


def print_ok(msg):
    """
    Prints the given message with an [OK] prefix.
    """
    print("\x1B[1;32m[OK]\x1B[0m    --> {}".format(msg))


