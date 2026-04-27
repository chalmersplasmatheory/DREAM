import h5py
import matplotlib.pyplot as plt
import numpy as np

from DREAM import DREAMOutput

do = DREAMOutput('output.h5')
time = do.grid.t[:]
ion_density = do.eqsys.n_i['D'][1][:,0]

plt.plot(time, ion_density)
plt.xlabel('Tid (s)')
plt.ylabel('Jon-densitet (n_i)')
plt.title('Utveckling av jon-densitet över tid')
plt.grid(True)
plt.show()


