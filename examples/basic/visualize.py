#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py')

from DREAM.DREAMOutput import DREAMOutput


do = None
if len(sys.argv) == 1:
    do = DREAMOutput('output.h5')
elif len(sys.argv) == 2:
    do = DREAMOutput(sys.argv[1])
else:
    print('ERROR: Invalid command line arguments. Expected at most one argument.')
    sys.exit(1)

# Declare unknowns (for convenience)
E_field = do.eqsys.E_field
n_cold  = do.eqsys.n_cold
n_hot   = do.eqsys.n_hot
n_re    = do.eqsys.n_re
T_cold  = do.eqsys.T_cold
n_tot   = do.eqsys.n_tot
n_i     = do.eqsys.n_i
f_hot   = do.eqsys.f_hot

#do.eqsys.E_field.plotRadialProfile(t=0)
#do.eqsys.E_field.plot()
#plt.show()

