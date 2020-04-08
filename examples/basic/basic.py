#!/usr/bin/env python3
#
# A very basic DREAM Python example. This script generates a basic
# DREAM input file which can be passed to 'dreami'.
#
# Run as
#
#   $ ./basic.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')
from DREAMSettings import DREAMSettings


ds = DREAMSettings()

ds.timestep.setTmax(1.0)
ds.timestep.setNt(20)

ds.save('dream_settings.h5')

