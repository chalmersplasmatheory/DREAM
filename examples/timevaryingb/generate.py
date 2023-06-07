#!/usr/bin/env python3
#
# This example illustrates how the time-varying B operator can be used.
#
# Run as
#
#   $ ./generate.py
#   $ ../../build/iface/dreami dream_settings.h5
#
# ###################################################################

from argparse import ArgumentParser
import numpy as np
import sys

sys.path.append('../../py/')

import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
from DREAM import DREAMSettings, runiface


def generate(dBdt=-1.0, B0=1.45, a=0.25):
    ds = DREAMSettings()
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_PARTIALLY_SCREENED

    # Physical parameters
    E = 0.1     # Electric field strength (V/m)
    n = 1e19    # Electron density (m^-3)
    T = 1000    # Temperature (eV)

    # Grid parameters
    pMax = 50   # maximum momentum in units of m_e*c
    Np   = 1000 # number of momentum grid points
    Nxi  = 50   # number of pitch grid points
    tMax = 1    # simulation time in seconds
    Nt   = 200  # number of time steps

    # Set E_field
    ds.eqsys.E_field.setPrescribedData(E)

    # Set temperature
    ds.eqsys.T_cold.setPrescribedData(T)

    # Set ions
    ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=n)

    # Disable avalanche generation
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)

    # Hot-tail grid settings
    ds.hottailgrid.setNxi(Nxi)
    ds.hottailgrid.setNp(Np)
    ds.hottailgrid.setPmax(pMax)

    # Set initial hot electron Maxwellian
    ds.eqsys.f_hot.setInitialProfiles(n0=n, T0=T)

    # Set boundary condition type at pMax
    ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_F_0) # F=0 outside the boundary
    ds.eqsys.f_hot.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)
    ds.eqsys.f_hot.setAdvectionInterpolationMethod(DistFunc.AD_INTERP_UPWIND)

    # Disable runaway grid
    ds.runawaygrid.setEnabled(False)

    # Set up radial grid
    ds.radialgrid.setB0(B0)
    ds.radialgrid.setMinorRadius(a)
    ds.radialgrid.setWallRadius(a*1.2)
    ds.radialgrid.setNr(1)

    if dBdt != 0:
        ds.radialgrid.setTimeVaryingB(dB0dt_B0=dBdt/B0)
        ds.eqsys.f_hot.setTimeVaryingB(True)

    # Set solver type
    ds.solver.setType(Solver.LINEAR_IMPLICIT) # semi-implicit time stepping
    ds.solver.preconditioner.setEnabled(False)

    # include otherquantities to save to output
    ds.other.include('fluid', 'hottail/Ap2')

    if dBdt != 0:
        ds.other.include('hottail/timevaryingb_Ap2')

    # Set time stepper
    ds.timestep.setTmax(tMax)
    ds.timestep.setNt(Nt)

    ds.output.setTiming(stdout=True, file=True)

    runiface(ds, f'output_dBdt_{dBdt:.3f}.h5', quiet=False)


def main():
    parser = ArgumentParser(description='Test time-varying B-field operator')
    parser.add_argument('-a', '--a', help="Tokamak minor radius.", type=float, nargs='?', default=0.25)
    parser.add_argument('-B', '--B0', help="On-axis magnetic field strength.", type=float, nargs='?', default=1.45)
    parser.add_argument('dBdt', help="Values of magnetic field time-rate-of-change to test.", nargs='*', type=float, default=[0.0, -1.0, -3.0, -6.0])
    args = parser.parse_args()

    for dBdt in args.dBdt:
        generate(dBdt, B0=args.B0, a=args.a)

    return 0


if __name__ == '__main__':
    sys.exit(main())


