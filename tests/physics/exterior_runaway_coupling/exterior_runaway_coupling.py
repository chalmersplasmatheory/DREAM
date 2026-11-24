"""Match exterior drag and ionization to the coefficients at the outer face."""

from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from scipy.constants import c, m_e

import dreamtests
from DREAM import DREAMSettings, runiface
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Solver as Solver


SEED = 1e17
DT = 1e-7
NXI = 8


def set_grid(grid, equation, pmin, pmax, np1=16):
    grid.setEnabled(True)
    grid.setNxi(NXI)
    grid.setNp(np1)
    grid.setPmin(pmin)
    grid.setPmax(pmax)
    p = pmin + (np.arange(np1) + .5)*(pmax-pmin)/np1
    xi = -1 + (np.arange(NXI) + .5)*2/NXI
    equation.setInitialValue(np.zeros((1,NXI,np1)), r=[.25], p=p, xi=xi)


def settings(layout, p, mode):
    ds = DREAMSettings()
    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)
    if layout in ('hot', 'both'):
        set_grid(ds.hottailgrid, ds.eqsys.f_hot, 0, p if layout == 'hot' else .4)
    if layout == 'both':
        set_grid(ds.runawaygrid, ds.eqsys.f_re, .4, p)
    elif layout == 'resolved':
        # The first cell center is exactly p, providing a kinetic v*sigma reference.
        set_grid(ds.runawaygrid, ds.eqsys.f_re, 0, 4*p, np1=2)

    if layout == 'hot':
        ds.eqsys.f_hot.setExternalKineticBoundary()
    elif layout == 'both':
        ds.eqsys.f_re.setExternalKineticBoundary()

    ds.radialgrid.setB0(5.)
    ds.radialgrid.setNr(1)
    ds.radialgrid.setMinorRadius(.5)
    ds.radialgrid.setWallRadius(.5)
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
    ds.collisions.collfreq_type = Collisions.COLLFREQ_TYPE_COMPLETELY_SCREENED
    ds.collisions.lnlambda = Collisions.LNLAMBDA_ENERGY_DEPENDENT
    ds.eqsys.E_field.setPrescribedData(0.)
    ds.eqsys.n_cold.setPrescribedData(5e19)
    ds.eqsys.T_cold.setType(Temperature.TYPE_SELFCONSISTENT)
    ds.eqsys.T_cold.setInitialProfile(1000.)
    ds.eqsys.T_cold.recombination = False
    ds.eqsys.n_i.addIon('H', Z=1, iontype=Ions.IONS_DYNAMIC, n=np.array([[1e18],[5e19]]), r=np.array([.25]), hydrogen=True)
    ds.eqsys.n_i.setIonization(mode)
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.timestep.setTmax(DT)
    ds.timestep.setNt(1)
    return ds


def simulate(ds, path):
    output = runiface(ds, outfile=str(path), quiet=True)
    try:
        np.testing.assert_allclose(output.grid.t, [0,DT], atol=0)
    finally:
        output.close()


def run(args):
    kinetic = Ions.IONIZATION_MODE_KINETIC
    with TemporaryDirectory() as directory:
        root = Path(directory)
        for p in (.8, 40.):
            ds = settings('resolved', p, kinetic)
            ds.other.include('runaway/kinioniz_vsigma')
            reference = root/'reference.h5'
            simulate(ds, reference)
            with h5py.File(reference) as f:
                vsigma = f['other/runaway/kinioniz_vsigma'][-1,0,0,0,0]
            assert vsigma > 0

            ds = settings('hot', p, kinetic)
            ds.eqsys.n_re.setInitialProfile(SEED)
            checkpoint = root/'checkpoint.h5'
            simulate(ds, checkpoint)
            with h5py.File(checkpoint, 'r+') as f:
                # Opposite exterior populations cancel in current and add in both rates.
                for pole in ('pos','neg'):
                    f[f'eqsys/n_external_{pole}'][-1,0] = SEED/2

            ignore = ['f_hot','j_re','j_tot','I_p']
            ds = settings('hot', p, Ions.IONIZATION_MODE_FLUID)
            ds.fromOutput(str(checkpoint), ignore=ignore)
            thermal = root/'thermal.h5'
            simulate(ds, thermal)
            with h5py.File(thermal) as f:
                neutral_thermal = f['eqsys/n_i'][-1,0,0]

            rates = []
            for layout, mode in (('hot', kinetic), ('both', kinetic), ('both', Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)):
                ds = settings(layout, p, mode)
                ds.fromOutput(str(checkpoint), ignore=ignore)
                grid = 'hottail' if layout == 'hot' else 'runaway'
                ds.other.include('fluid/Tcold_nre_coll', 'fluid/Tcold_radiation',
                                 'fluid/Tcold_binding_energy',
                                 'fluid/reioniz_rate', 'fluid/reioniz_vsigma', f'{grid}/nu_s_f1')
                result = root/'result.h5'
                simulate(ds, result)
                with h5py.File(result) as f:
                    nu_s = f[f'other/{grid}/nu_s_f1'][-1,0,0,-1]
                    heat = -f['other/fluid/Tcold_nre_coll'][-1,0]
                    rad = f['other/fluid/Tcold_radiation'][-1,0]
                    binding = f['other/fluid/Tcold_binding_energy'][-1,0]
                    w = f['eqsys/W_cold'][:,0]
                    expected_heat = m_e*c*p * (c*p/np.hypot(1,p)) * nu_s * SEED
                    # Allow for differences between SciPy's and DREAM's electron mass constants.
                    np.testing.assert_allclose(heat, expected_heat, rtol=2e-7)
                    np.testing.assert_allclose((w[-1]-w[0])/DT + rad + binding, heat, rtol=1e-7)
                    coeff = f['other/fluid/reioniz_vsigma'][-1,0,0]
                    np.testing.assert_allclose(coeff, vsigma*SEED, rtol=1e-10)
                    rate = f['other/fluid/reioniz_rate'][-1,:,0]
                    neutral = f['eqsys/n_i'][-1,0,0]
                    np.testing.assert_allclose(rate, np.array([-1,1]) * coeff*neutral, rtol=1e-10)
                    # The extra neutral loss responds to the impact term as well as
                    # the consequent change in thermal rates over this step.
                    np.testing.assert_allclose((neutral-neutral_thermal)/DT, rate[0], rtol=.15)
                    np.testing.assert_allclose(f['eqsys/j_re'][-1,0], 0, atol=1e-8)
                    rates.append((heat, neutral))
            np.testing.assert_allclose(rates[1:], [rates[0], rates[0]], rtol=1e-10)
            if args['verbose']:
                dreamtests.print_ok(f'p_max={p:g}: heating {rates[0][0]:.6g} W/m^3, v*sigma {vsigma:.6g} m^3/s.')
    dreamtests.print_ok('Exterior heating and ionization match the kinetic outer-face coefficients.')
    return True
