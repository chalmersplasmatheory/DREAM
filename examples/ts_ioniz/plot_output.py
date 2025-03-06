import h5py
import matplotlib.pyplot as plt
import numpy as np

with h5py.File('output.h5', 'r') as f:
    time = np.array(f['grid/t'])
    ion_density = np.array(f['eqsys/n_i'])

    plt.plot(time, ion_density)
    plt.xlabel('Tid (s)')
    plt.ylabel('Jon-densitet (n_i)')
    plt.title('Utveckling av jon-densitet över tid')
    plt.grid(True)
    plt.show()


