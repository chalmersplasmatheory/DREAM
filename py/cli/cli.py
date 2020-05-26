#!/usr/bin/python3 -i

import matplotlib.pyplot as plt
import numpy as np
import pathlib
import sys

# TODO Set this up properly
sys.path.append(pathlib.Path(__file__).parent.parent.absolute())

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput


do = None
if len(sys.argv) == 1:
    do = DREAMOutput('output.h5')
elif len(sys.argv) == 2:
    do = DREAMOutput(sys.argv[1])
else:
    print('ERROR: Invalid command line arguments. Expected at most one argument (name of file).')
    sys.exit(1)

setup_interactive(do, glob=globals())

