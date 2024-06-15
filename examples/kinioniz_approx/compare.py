#!/usr/bin/env python3
#
# Compares the fluid RE impact ionization mode with the fully kinetic one.
#
# Run as
#
#   $ ./compare.py
#
# And then run plot script as
#
#   $ ./plot.py
#
# ###################################################################

import sys, os
import numpy as np

import DREAM
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.RunawayElectronDistribution as REDist

sys.path.append("../../py/DREAM/Formulas/")
from Distributions import getAvalancheDistribution

ARGON_DENSITY = 1e22

def run(ds, tmax, outfile="output.h5"):
    if os.path.exists(outfile):
        return DREAM.DREAMOutput(outfile)
    ds.timestep.setIonization(dt0=1e-7, dtmax=1e-6, tmax=tmax)
    return DREAM.runiface(ds, outfile=outfile)

def generate_base(nre, temperature, electric_field):

    ds = DREAM.DREAMSettings()

    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(1)
    ds.radialgrid.setWallRadius(1.1)


    if isinstance(nre, (int, float)):
        ds.eqsys.n_re.setInitialProfile(nre)
        ds.radialgrid.setNr(1)
    elif isinstance(nre, (list, np.ndarray)):
        ds.eqsys.n_re.setInitialProfile(nre, np.linspace(0, 1, len(nre)))
        ds.radialgrid.setNr(len(nre))

    ds.eqsys.E_field.setPrescribedData(electric_field)
    ds.eqsys.T_cold.setPrescribedData(temperature)
    ds.eqsys.n_i.addIon(name='Ar', Z=18, n=ARGON_DENSITY, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED)

    return ds


def generate_fluid(nre, temperature, electric_field):
    ds = generate_base(nre, temperature, electric_field)
    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_FLUID_RE)
    ds.other.include(["fluid/Zeff", "fluid/kinioniz_vsigma"])
    return ds


def generate_kinetic(nre, temperature, electric_field):
    ds = generate_base(nre, temperature, electric_field)
    assert ds.radialgrid.nr == 1
    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC)

    Np = 100
    Nxi = 40
    pMin, pMax = .1, 100

    ds.runawaygrid.setEnabled(True)
    ds.runawaygrid.setNp(Np)
    ds.runawaygrid.setNxi(Nxi)
    ds.runawaygrid.setPmin(pMin)
    ds.runawaygrid.setPmax(pMax)

    f   = np.zeros((1, 1, Nxi, Np))
    p  = np.linspace(pMin, pMax, Np)
    xi = np.linspace(-1, 1, Nxi)

    f[0,0,:] = getAvalancheDistribution(p=p, xi=xi, E=electric_field, Z=18, nre=nre)

    ds.eqsys.f_re.prescribe(f=f, t=[0], r=[0], xi=xi, p=p)

    ds.eqsys.f_re.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)
    ds.other.include(["fluid/Zeff", "runaway/kinioniz_vsigma"])
    return ds


if __name__ == '__main__':


    LOG_NRE_MIN = 15
    LOG_NRE_MAX = 19
    NSCAN_FLUID = 100
    NSCAN_KINETIC = 6

    import argparse

    parser = argparse.ArgumentParser(description="Run scan in 'n_re' to compare RE impact ionization between two models.")
    parser.add_argument("-T", "--temperature", help="Temperature of both electrons (and ions) [eV].", dest="temperature", action="store", type=float)
    parser.add_argument("-E", "--electric_field", help="Electric field (only relevant for the analytical avalanche runaway distribution) [m^-3]", dest="electric_field", action="store", type=float)
    parser.add_argument("-d", "--save_dir", help="Name of directory to save outputs to.", dest="save_dir", action="store", type=str)
    parser.set_defaults(temperature=1, electric_field=1, save_dir="outputs")

    args = parser.parse_args()

    from pathlib import Path
    Path(args.save_dir).mkdir(parents=True, exist_ok=True)

    ds = generate_fluid(np.logspace(LOG_NRE_MIN, LOG_NRE_MAX, NSCAN_FLUID), temperature=args.temperature, electric_field=args.electric_field)
    do = run(ds, 1e-3, outfile=f"{args.save_dir}/fluid.h5")

    for i, nre in enumerate(np.logspace(LOG_NRE_MIN, LOG_NRE_MAX, NSCAN_KINETIC)):
        ds = generate_kinetic(nre, temperature=args.temperature, electric_field=args.electric_field)
        do = run(ds, 1e-3, outfile=f"{args.save_dir}/kinetic_{i}.h5")
