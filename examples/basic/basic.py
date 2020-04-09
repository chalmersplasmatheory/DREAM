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

import Solver


ds = DREAMSettings()

# Set equations
times  = [0]
radius = [0, 1]
density = 1e19 * np.ones((len(times), len(radius)))
ds.equationsystem.n_cold.setPrescribedData(density=density, times=times, radius=radius)

# Hot-tail grid settings
ds.hottailgrid.setNp(100)
ds.hottailgrid.setPmax(5)

# Disable runaway grid
ds.runawaygrid.setEnabled(False)

# Set up radial grid
ds.radialgrid.setB0(5)
ds.radialgrid.setMinorRadius(0.22)
ds.radialgrid.setNr(10)

# Use the linear solver
ds.solver.setType(Solver.TYPE_LINEAR_IMPLICIT)

# Set time stepper
ds.timestep.setTmax(1.0)
ds.timestep.setNt(20)

# Save settings to HDF5 file
ds.save('dream_settings.h5')

