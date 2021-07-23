#!/usr/bin/python3 -i

import sys

sys.path.append('../../py/')
sys.path.append('../../')
sys.path.append('../../build/dreampyface/cxx/')

import dreampyface
from DREAM import DREAMSettings, DREAMOutput


def callback(ptr):
    sim = dreampyface.Simulation(ptr)
    print('TIME = {}'.format(sim.getCurrentTime()))

ds = DREAMSettings()

n = 1e19
T = 1e3

ds.eqsys.E_field.setPrescribedData(1)
ds.eqsys.n_i.addIon(name='D', Z=1, n=n)
ds.eqsys.T_cold.setPrescribedData(T)
ds.eqsys.f_hot.setInitialProfiles(rn0=0, n0=n, rT0=0, T0=T)

ds.hottailgrid.setNxi(20)
ds.hottailgrid.setNp(500)
ds.hottailgrid.setPmax(5)

ds.runawaygrid.setEnabled(False)

ds.radialgrid.setB0(1)
ds.radialgrid.setMinorRadius(0.1)
ds.radialgrid.setNr(1)
ds.radialgrid.setWallRadius(0.12)

ds.timestep.setTmax(1e-3)
ds.timestep.setNt(20)

ds.output.setFilename('output.h5')

dreampyface.register_callback_timestep_finished(callback)

s = dreampyface.setup_simulation(ds)

unknowns = s.getUnknowns()
for uq, info in unknowns.items():
    print('{:12s}  {:8d}  {}'.format(uq, info['nelements'], info['description']))

#do = s.run()

