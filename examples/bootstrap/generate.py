#!/usr/bin/env python3
#
# A simple example of how to set up a DREAM simulation with bootstrap current.
#
import numpy as np

from DREAM.DREAMSettings import DREAMSettings

import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.OhmicCurrent as OhmicCurrent
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.BootstrapCurrent as BootstrapCurrent
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature

import ASDEXU

NT = 1
TMAX = 1e-3

NR = 100

INCLUDE_BOOTSTRAP_IN_INITIAL_CURRENT = False
NONLINEAR_SOLVER = True

def getProfiles():
    """
    Returns profiles for electron temperature, electron density and ion temperature.
    These profiles are fitted from AUGD/IDF data for ASDEX Upgrade discharge #32305 @ t = 1.8103 s.
    """
    def p(p0, coeffs):
        return lambda r: p0*np.poly1d(np.flip(coeffs))(r**2)

    te_coeffs = [0.994, -0.922, 19.837, -211.339, 900.765, -1978.44, 2360.032, -1449.313, 358.447]
    ne_coeffs = [0.994, 0.941, -24.927, 206.672, -884.368, 2065.323, -2670.298, 1797.844, -491.933]
    ti_coeffs = [0.997, -5.179, 26.581, -96.837, 240.105, -378.705, 349.822, -167.213, 30.567]

    return p(4.235e3, te_coeffs), p(5.466e19, ne_coeffs), p(4.019e3, ti_coeffs)

def generate():
    """
    Sets up a simple DREAM simulation with bootstrap current included and returns the corresponding DREAMSettings object.
    """

    ds = DREAMSettings()

    ASDEXU.setMagneticField(ds, nr=NR)
    a = ASDEXU.a

    r = np.linspace(0, 1, NR)

    te, ne, ti = getProfiles()


    # Set electric field
    ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
    ds.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_PRESCRIBED, V_loop_wall_R0=0, R0=ds.radialgrid.R0)


    # Set Ohmic current
    johm = 1 - .9 * r**2
    ds.eqsys.j_ohm.setInitialProfile(johm, radius=r*a, Ip0=ASDEXU.Ip)
    ds.eqsys.j_ohm.setConductivityMode(OhmicCurrent.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)

    # Include bootstrap current
    ds.eqsys.j_bs.setMode(BootstrapCurrent.BOOTSTRAP_MODE_ENABLED)
    if INCLUDE_BOOTSTRAP_IN_INITIAL_CURRENT:
        ds.eqsys.j_bs.setInitMode(BootstrapCurrent.BOOTSTRAP_INIT_MODE_TOTAL)
    else:
        ds.eqsys.j_bs.setInitMode(BootstrapCurrent.BOOTSTRAP_INIT_MODE_OHMIC) # default init mode

    # Set temperature
    ds.eqsys.T_cold.setPrescribedData(te(r), radius=r*a)
    ds.eqsys.T_cold.setType(Temperature.TYPE_PRESCRIBED)


    # Add ions
    ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=.8*ne(r), r=r*a, T=ti(r))
    ds.eqsys.n_i.addIon(name='B', Z=5, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=.04*ne(r), r=r*a, T=ti(r))

    # Disable kinetic grids
    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    # Select solver
    if NONLINEAR_SOLVER:
        ds.solver.setType(Solver.NONLINEAR)
    else:
        ds.solver.setType(Solver.LINEAR_IMPLICIT)

    # Set time stepper
    ds.timestep.setTmax(TMAX)
    ds.timestep.setNt(NT)

    return ds


if __name__ == '__main__':

    ds = generate()
    ds.save("dream_settings.h5")
