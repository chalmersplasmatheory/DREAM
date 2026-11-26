"""Check runaway generation and pitch reflection under E -> -E."""

import numpy as np

import dreamtests

from DREAM import DREAMSettings, runiface
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.DistributionFunction as Distribution
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver


E_FIELD = 3.0  # V/m
T_MAX = 4e-4  # s
N_STEPS = 16
P_MIN = 0.16
SEEDS = ('gammaDreicer', 'gammaTritium', 'gammaCompton')
SOURCES = ('avalanche', 'Dreicer', 'tritium', 'Compton')


def settings(toroidal, direction):
    ds = DREAMSettings()
    ds.hottailgrid.setEnabled(False)
    grid = ds.runawaygrid
    # Resolve tritium births below p ~ 0.27 while retaining a substantial
    # population of independently passing Rosenbluth–Putvinski knock-ons.
    grid.setPmin(P_MIN)
    grid.setPmax(3)
    grid.setNp(56)
    grid.setBiuniformGrid(psep=0.3, npsep=20)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(0.4)
    ds.radialgrid.setWallRadius(0.4)
    ds.radialgrid.setNr(1)
    if toroidal:
        ds.radialgrid.setMajorRadius(2)
        ds.radialgrid.setShaping(psi=[0, 0.02], rpsi=[0, 0.4], GOverR0=1)
        grid.setTrappedPassingBoundaryLayerGrid(
            dxiMax=0.1, NxiPass=12, NxiTrap=6, boundaryLayerWidth=0.004
        )
    else:
        grid.setNxi(32)

    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
    ds.eqsys.E_field.setPrescribedData(direction*E_FIELD)
    ds.eqsys.T_cold.setPrescribedData(45)
    ds.eqsys.n_cold.setPrescribedData(2e19)
    for species in ('D', 'T'):
        ds.eqsys.n_i.addIon(species, Z=1, n=1e19, tritium=(species == 'T'),
                           iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED)

    re = ds.eqsys.n_re
    re.setInitialProfile(0)
    re.setNegativeRunaways()
    re.setAvalanche(Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=P_MIN)
    re.setTritium(Runaways.TRITIUM_MODE_KINETIC)
    re.setCompton(Runaways.COMPTON_RATE_ITER_DMS_KINETIC)
    re.setDreicer(Runaways.DREICER_RATE_CONNOR_HASTIE)
    ds.eqsys.f_re.setAdvectionInterpolationMethod(ad_int=Distribution.AD_INTERP_UPWIND)
    ds.eqsys.f_re.setBoundaryCondition(Distribution.BC_F_0)

    # The two-way boundary is transport-only. This fixture exercises the
    # source operators on a resolved runaway grid.
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.timestep.setTmax(T_MAX)
    ds.timestep.setNt(N_STEPS)
    ds.other.include(*(f'fluid/{seed}' for seed in SEEDS), 'fluid/runawayRate',
                     'runaway/S_ava', 'runaway/S_compton', 'runaway/S_tritium')
    return ds


def read_case(toroidal, direction):
    out = runiface(settings(toroidal, direction), quiet=True)
    try:
        t = np.asarray(out.grid.t)
        np.testing.assert_allclose(t, np.linspace(0, T_MAX, N_STEPS+1), rtol=1e-12, atol=0)
        np.testing.assert_allclose(out.eqsys.E_field[:], direction*E_FIELD, rtol=1e-12)
        f = out.eqsys.f_re
        grid = f.momentumgrid
        weights = grid.Vprime_VpVol[0]*grid.DP2[0]*grid.DP1[0]
        faces = np.asarray(grid.p2_f)
        np.testing.assert_allclose(faces, -faces[::-1], rtol=0, atol=1e-12)
        trapped = float(out.grid.xi0TrappedBoundary[0])

        fluid = out.other.fluid
        kinetic = out.other.runaway
        seeds = {name: fluid[name][:, 0] for name in SEEDS}

        # S_ava evaluates the f_re, n_re block. At the first momentum
        # cell the RP source lies well away from the pitch poles, so the
        # two corner cells measure the Dreicer deposition directly.
        ava = -kinetic.S_ava[:, 0, :, :]*weights
        dreicer = np.zeros_like(ava)
        dreicer[:, [0, -1], 0] = ava[:, [0, -1], 0]
        ava[:, [0, -1], 0] = 0
        pole = -1 if direction > 0 else 0
        opposite = 0 if direction > 0 else -1
        np.testing.assert_allclose(dreicer[:, pole, 0], seeds['gammaDreicer'], rtol=1e-10)
        np.testing.assert_allclose(dreicer[:, opposite, 0], 0,
                                   atol=1e-10*np.max(seeds['gammaDreicer']))
        sources = {
            'avalanche': ava,
            'Dreicer': dreicer,
            'tritium': kinetic.S_tritium[:, 0, :, :]*weights,
            'Compton': kinetic.S_compton[:, 0, :, :]*weights,
        }
        rates = {
            **seeds,
            'dn_re/dt': fluid.runawayRate[:, 0],
            'n_re': out.eqsys.n_re[1:, 0],
        }
        return faces, trapped, weights, sources, rates, f[:, 0, :, :], out.eqsys.j_re[:, 0]
    finally:
        out.close()


def relative_error(a, b):
    """Maximum difference in a history, scaled by the peak magnitude."""
    return np.max(np.abs(a-b))/max(np.max(np.abs(a)), np.max(np.abs(b)))


def reflected_error(plus, minus, passing):
    """Peak L1 mismatch of cell populations in the independent passing region."""
    a = plus[:, passing, :]
    b = minus[:, ::-1, :][:, passing, :]
    norm = max(np.sum(np.abs(x), axis=(1, 2)).max() for x in (a, b))
    return np.sum(np.abs(a-b), axis=(1, 2)).max()/norm


def check(toroidal, verbose):
    faces, trapped, weights, positive, rates_pos, f_pos, j_pos = read_case(toroidal, +1)
    faces_neg, trapped_neg, weights_neg, negative, rates_neg, f_neg, j_neg = read_case(toroidal, -1)
    np.testing.assert_array_equal(faces, faces_neg)
    np.testing.assert_allclose(trapped, trapped_neg)
    np.testing.assert_allclose(weights, weights_neg, rtol=1e-12)
    if (trapped > 0) != toroidal:
        raise AssertionError("The toroidal fixture must have a trapped region.")

    # Toroidal trapped rows share orbit populations. The negative cell
    # straddling -xi_T can have zero measure while its positive counterpart
    # retains weight. Use all native weights in scalar moments, and use
    # only fully passing cells for cell-by-cell reflection.
    passing = (faces[1:] <= -trapped) | (faces[:-1] >= trapped)
    np.testing.assert_array_equal(passing, passing[::-1])
    if not np.any(passing):
        raise AssertionError("The fixture must contain independently resolved passing cells.")

    # The cylindrical passing distribution has a reflected partner in each
    # pitch cell. Toroidal reflection and current allow the measured
    # sensitivity to the separatrix discretization.
    scalar_tol = 1e-5
    source_pitch_tol = 1e-2 if toroidal else 2e-4
    distribution_tol = 1e-3 if toroidal else 1e-4
    current_tol = 1e-3 if toroidal else 1e-4
    results = {}

    for name in rates_pos:
        a, b = rates_pos[name], rates_neg[name]
        if any(not np.isfinite(x).all() or np.max(x) <= 0 for x in (a, b)):
            raise AssertionError(f"{name} must be finite and active for both field signs.")
        results[name] = (relative_error(a, b), scalar_tol)

    for name in SOURCES:
        a, b = positive[name], negative[name]
        moments = (a.sum(axis=(1, 2)), b.sum(axis=(1, 2)))
        if any(not np.isfinite(x).all() or np.max(x) <= 0 for x in moments):
            raise AssertionError(f"{name} must have finite, positive integrated production for both field signs.")
        results[f'{name} kinetic moment'] = (relative_error(*moments), scalar_tol)
        results[f'{name} passing reflection'] = (reflected_error(a, b, passing), source_pitch_tol)

    a, b = f_pos*weights, f_neg*weights_neg
    results['resolved density'] = (relative_error(a.sum(axis=(1, 2)), b.sum(axis=(1, 2))), scalar_tol)
    results['f_re passing reflection'] = (reflected_error(a, b, passing), distribution_tol)
    if not (j_pos[-1] > 0 and j_neg[-1] < 0):
        raise AssertionError("The runaway current must follow the field direction.")
    results['j_re sign reversal'] = (relative_error(j_pos, -j_neg), current_tol)

    label = 'Toroidal' if toroidal else 'Cylindrical'
    if verbose and toroidal:
        j = np.searchsorted(faces, -trapped)-1
        print(f"    xi_T = {trapped:.6f}, separatrix cell measures "
              f"(negative, positive) = ({weights[j].sum():.6g}, {weights[-1-j].sum():.6g})")

    success = True
    for name, (error, tolerance) in results.items():
        message = f'{label} {name}, error {error:.3g}, tolerance {tolerance:.3g}'
        if np.isfinite(error) and error <= tolerance:
            if verbose:
                dreamtests.print_ok(message)
        else:
            dreamtests.print_error(message)
            success = False
    if success:
        dreamtests.print_ok(f'{label} field-reversal symmetry with all four sources enabled.')
    return success


def run(args):
    cylinder = check(toroidal=False, verbose=args['verbose'])
    torus = check(toroidal=True, verbose=args['verbose'])
    return cylinder and torus
