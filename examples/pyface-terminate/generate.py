#!/usr/bin/env python3

from DREAM import DREAMSettings, runiface
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../build/dreampyface/cxx/')
sys.path.append('../../')
import dreampyface

import DREAM.Settings.Equations.ElectricField as EField
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as RE
import DREAM.Settings.Solver as Solver


def terminate_ioniz(sim):
    """
    Function which determines when to stop the ionization simulation.
    """
    # Fractional change in ncold below which ionization should stop
    THRESHOLD = 1e-4

    ncold = sim.unknowns.getData('n_cold')
    ntot  = sim.unknowns.getData('n_tot')

    if ncold['x'].shape[0] < 2:
        return False

    dnc = ncold['x'][-1,:] - ncold['x'][-2,:]

    mx = np.amax(dnc/ntot['x'][-1,:])
    #print(f'Max change: {mx}')

    return mx < THRESHOLD


def baseline():
    """
    Generate a simulation baseline settings object.
    """
    ds = DREAMSettings()

    a = 0.5
    b = 0.55
    Ip = 800e3

    ds.radialgrid.setB0(3.1)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(b)
    ds.radialgrid.setNr(30)

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    ds.eqsys.T_cold.setPrescribedData(5.8e3)
    ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=2.6e19)
    ds.eqsys.E_field.setPrescribedData(0)

    ds.solver.setType(Solver.LINEAR_IMPLICIT)
    ds.solver.setLinearSolver(linsolv=Solver.LINEAR_SOLVER_MKL)

    ds.timestep.setTmax(1)
    ds.timestep.setNt(1)

    ds.other.include('fluid/conductivity')

    do = runiface(ds)

    ###########################
    # Calculate electric field
    ###########################
    sigma = do.other.fluid.conductivity[-1,:]

    jprof = np.ones((do.grid.r[:].size,))
    j0 = Ip / (2*np.pi*do.grid.integrate(jprof * do.grid.r[:]))
    E = j0*jprof / sigma

    ds.eqsys.E_field.setType(EField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setInitialProfile(E, radius=do.grid.r[:])
    ds.eqsys.E_field.setBoundaryCondition(EField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0)

    ds.eqsys.n_i.addIon('Ar', Z=18, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e20)

    ds.timestep.setTerminationFunction(terminate_ioniz)
    ds.timestep.setNt(2000)
    ds.timestep.setTmax(1e-4)

    do.close()

    return ds


def main():
    ds = baseline()

    #runiface(ds, 'output_ioniz.h5')
    s = dreampyface.setup_simulation(ds)
    do = s.run()

    dnc = np.diff(do.eqsys.n_cold[:,0])
    plt.plot(do.grid.t[1:], dnc/do.eqsys.n_tot[-1,0])
    plt.show()

    return 0


if __name__ == '__main__':
    sys.exit(main())


