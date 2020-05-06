#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py')

from DREAMOutput import DREAMOutput


do = DREAMOutput('output.h5')

do.eqsys.E_field.plotRadialProfile(t=0)
plt.show()
