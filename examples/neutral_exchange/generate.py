#!/usr/bin/env python3
"""Isolate neutral D--Ne heat exchange at fixed particle densities."""

from pathlib import Path
import sys

import numpy as np
from scipy.constants import e

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / 'py'))

from DREAM import DREAMSettings, DREAMIO, runiface
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


def run(order, label):
    ds = DREAMSettings()
    ds.eqsys.E_field.setPrescribedData(0)
    ds.eqsys.T_cold.setPrescribedData(10)

    # Two fixed neutral populations. A fixed D+ background supplies electrons
    # for the surrounding DREAM fluid equations; it does not exchange heat
    # with W_n in the present implementation.
    neutral_density = {'D': 1e19, 'Ne': 1e19}
    temperature = {'D': 2.0, 'Ne': 8.0}  # eV; equilibrium is 6 eV
    charge = {'D': 1, 'Ne': 10}
    radius = np.array([0.0, 0.42])

    for name in order:
        density = np.zeros((charge[name] + 1, 1, radius.size))
        density[0, 0, :] = neutral_density[name]
        if name == 'D':
            density[1, 0, :] = 1e18
        ds.eqsys.n_i.addIon(
            name=name, Z=charge[name], iontype=Ions.IONS_PRESCRIBED,
            n=density, r=radius, t=np.array([0.0]), T=10
        )

    # Prescribed densities disable particle evolution. No reactions, particle
    # sources, or transport are enabled.
    ds.eqsys.n_i.setNeutralTemperatureEnabled(True)
    energy = np.array([
        [1.5 * e * neutral_density[name] * temperature[name]]
        for name in order
    ])
    ds.eqsys.W_n.setInitialProfile(energy=energy, radius=np.array([0.0]))

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)
    ds.radialgrid.setB0(5)
    ds.radialgrid.setMinorRadius(0.42)
    ds.radialgrid.setWallRadius(0.44)
    ds.radialgrid.setNr(1)

    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setVerbose(False)
    ds.timestep.setTmax(0.05)
    ds.timestep.setNt(2000)
    ds.output.setFilename(str(HERE / f'output_{label}.h5'))
    settings = ds.todict()
    # This build's HDF5 reader crashes on the empty reaction-name string.
    # Reactions are disabled, so omit the names and use the C++ default.
    del settings['eqsys']['n_i']['reactions']['names']
    settings_file = str(HERE / f'settings_{label}.h5')
    DREAMIO.SaveDictAsHDF5(settings_file, settings)
    output = runiface(settings_file, str(HERE / f'output_{label}.h5'), quiet=False)
    output.close()
    print(f'Finished {label}: species order {order}')


if __name__ == '__main__':
    run(('D', 'Ne'), 'forward')
    run(('Ne', 'D'), 'reverse')
    print('Run plot.py to plot and check both outputs.')
