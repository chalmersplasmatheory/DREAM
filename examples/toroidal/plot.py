import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.DREAMOutput import DREAMOutput
plt.rcParams.update({'font.size': 14})

do = DREAMOutput('output.h5')

t = do.grid.t[1:]
n_t = do.other.fluid.runawayRate[:,0]
plt.plot(t, n_t)
plt.savefig('testfigure.pdf')
