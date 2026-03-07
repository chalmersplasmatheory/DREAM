#!/usr/bin/env python3
#
# Plot solutions
# ######

from pathlib import Path

import matplotlib.pyplot as plt
import sys

from generate import RESTART_OUTPUT_NAME

sys.path.append('../../py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput

output = Path(__file__).parent / RESTART_OUTPUT_NAME
if not output.exists():
    raise Exception("Example has not been run yet. First run: generate.py; dreami init_settings.h5; dreami restart_settings.h5")

do = DREAMOutput(output)

all_radii = list(range(do.grid.r.size))
plt.figure()
do.eqsys.E_field.plot(r=all_radii)

plt.figure()
do.eqsys.T_cold.plot(r=all_radii)

plt.figure()
do.eqsys.n_re.plot(r=all_radii)

plt.figure()
do.eqsys.f_hot.plot(t=range(do.grid.t.size))

plt.show()
