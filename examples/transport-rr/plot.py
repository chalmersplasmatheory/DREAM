#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py')

from DREAM.DREAMOutput import DREAMOutput


do = DREAMOutput('output.h5')

do.eqsys.T_cold.plot(t=[0,5,10,15,-1], show=False)
plt.show()

