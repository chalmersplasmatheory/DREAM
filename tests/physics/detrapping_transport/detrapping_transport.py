"""Settings, executable, diagnostic and restart tests for finite-detrapping RR."""
from pathlib import Path
import h5py
import subprocess
import tempfile

import numpy as np
import dreamtests
from DREAM import DREAMOutput, DREAMSettings, runiface
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.HotElectronDistribution as FHot
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.RadialGrid as RadialGrid
import DREAM.Settings.Solver as Solver
import DREAM.Settings.TransportSettings as Transport


def settings(model='eigenmode'):
    ds = DREAMSettings()
    ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_SUPERTHERMAL
    ds.radialgrid.setType(RadialGrid.TYPE_ANALYTIC_TOROIDAL)
    ds.radialgrid.setMajorRadius(3)
    ds.radialgrid.setMinorRadius(0.5)
    ds.radialgrid.setWallRadius(0.5)
    ds.radialgrid.setNr(3)
    r = np.linspace(0, 0.5, 20)
    ds.radialgrid.setShaping(GOverR0=5, psi=0.1*r*r, rpsi=r)
    ds.eqsys.E_field.setPrescribedData(0.01)
    ds.eqsys.T_cold.setPrescribedData(100)
    ds.eqsys.n_cold.setPrescribedData(5e19)
    ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED_FULLY_IONIZED, n=5e19)
    ds.hottailgrid.setNxi(1)
    ds.hottailgrid.setNp(24)
    ds.hottailgrid.setPmax(1)
    ds.runawaygrid.setEnabled(False)
    ds.eqsys.n_re.setDreicer(False)
    ds.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
    ds.eqsys.n_re.setCompton(False)
    ds.eqsys.n_re.setTritium(False)
    ds.eqsys.n_re.setHottail(False)
    ds.eqsys.f_hot.setInitialProfiles(n0=0.99*5e19, T0=100)
    ds.eqsys.f_hot.setParticleSource(FHot.PARTICLE_SOURCE_EXPLICIT)
    ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)
    # Time dependence exercises the current branch's prescribed interpolator.
    dBB = np.array([[0.003], [0.006]])
    if model == 'rr':
        ds.eqsys.f_hot.transport.setMagneticPerturbation(dBB, t=[0, 2e-7], r=[0])
    else:
        ds.eqsys.f_hot.transport.setDetrappingRechesterRosenbluth(
            dBB, t=[0, 2e-7], r=[0]
        )
    ds.solver.setType(Solver.NONLINEAR)
    ds.solver.setLinearSolver(Solver.LINEAR_SOLVER_LU)
    ds.timestep.setTmax(2e-7)
    ds.timestep.setNt(2)
    ds.other.include('hottail/Drr', 'hottail/Drr_detrapping', 'hottail/kappa_detrapping',
                     'hottail/nu_D_f2', 'scalar/radialloss_f_hot')
    return ds


def check_settings():
    ts = Transport.TransportSettings(kinetic=True)
    ts.setDetrappingRechesterRosenbluth(0.003)
    saved = ts.todict()
    assert saved['detrapping_mode'] == 3
    assert 'detrapping_scale' not in saved
    for mode in ('eigenmode', 'kappa', 'collisional_eigenmode', 3):
        ts.setDetrappingRechesterRosenbluth(0.003, mode=mode, detrapping_scale=1)
        ts.verifySettings()
        restored = Transport.TransportSettings(kinetic=True)
        restored.fromdict(ts.todict())
        restored.verifySettings()
        assert restored.type == Transport.TRANSPORT_RECHESTER_ROSENBLUTH_DETRAPPING
        np.testing.assert_array_equal(restored.dBB, ts.dBB)

    for mode, scale, kinetic in [
        ('local', 1, True), ('generalized', 1, True), (1, 1, True), (2, 1, True),
        (3.5, 1, True), ('eigenmode', 2, True), ('eigenmode', 0, True),
        ('eigenmode', float('nan'), True), ('eigenmode', float('inf'), True),
        ('eigenmode', 1, False)
    ]:
        try:
            ts = Transport.TransportSettings(kinetic=kinetic)
            ts.setDetrappingRechesterRosenbluth(0.003, mode=mode, detrapping_scale=scale)
            ts.verifySettings()
        except Transport.TransportException:
            continue
        raise AssertionError(f'Accepted invalid settings: {mode}, {scale}, {kinetic}')


def check_rejected_file(executable, path, message):
    rejected = subprocess.run([str(executable), str(path)], capture_output=True, text=True)
    assert rejected.returncode != 0, f'Accepted {path.name}'
    assert message in rejected.stdout + rejected.stderr, rejected.stdout + rejected.stderr


def run(args):
    quiet = not args.get('verbose', False)
    try:
        check_settings()
        with tempfile.TemporaryDirectory(prefix='dream-detrapping-') as directory:
            directory = Path(directory)
            executable = Path(__file__).resolve().parents[3] / 'build/iface/dreami'
            outputs = {}
            for model in ('rr', 'eigenmode'):
                ds = settings(model)
                out = directory / f'{model}.h5'
                ds.output.setFilename(str(out))
                path = directory / f'{model}_settings.h5'
                ds.save(str(path))
                if model == 'eigenmode':
                    # Check compatibility with existing unit-scale eigenmode files.
                    with h5py.File(path, 'r+') as raw:
                        raw['eqsys/f_hot/transport/detrapping_scale'] = np.array([1.0])
                loaded = DREAMSettings(str(path))
                assert loaded.eqsys.f_hot.transport.type == ds.eqsys.f_hot.transport.type
                outputs[model] = runiface(loaded, str(out), quiet=quiet)

            rr = outputs['rr'].other.hottail.Drr[:]
            eigen = outputs['eigenmode'].other.hottail
            D = eigen.Drr[:]
            cap = eigen.Drr_detrapping[:]
            kappa = eigen.kappa_detrapping[:]
            assert np.isfinite(D).all() and (D >= 0).all()
            assert (D <= rr*(1+1e-10)).all() and (D[:, 1:] < rr[:, 1:]).any()
            np.testing.assert_array_equal(D[:, 0], rr[:, 0])
            np.testing.assert_array_equal(cap[:, 0], 0)
            assert (cap[:, 1:] > 0).all()
            np.testing.assert_allclose(D[:, 1:], (rr*cap/(rr+cap))[:, 1:], rtol=1e-10)
            assert kappa.shape == (2, 4) and (kappa[:, 1:] > 0).all()
            np.testing.assert_array_equal(kappa[:, 0], 0)
            np.testing.assert_array_equal(kappa[0], kappa[1])
            np.testing.assert_array_equal(outputs['rr'].other.hottail.Drr_detrapping[:], 0)
            np.testing.assert_array_equal(outputs['rr'].other.hottail.kappa_detrapping[:], 0)
            np.testing.assert_allclose(rr[1], (16/9)*rr[0], rtol=1e-12)
            assert np.isfinite(outputs['eigenmode'].other.scalar.radialloss_f_hot[:]).all()

            # Restart with zero drive at the first step and the original final drive at the second.
            ds = settings()
            ds.fromOutput(str(directory / 'eigenmode.h5'))
            ds.eqsys.f_hot.transport.setDetrappingRechesterRosenbluth(
                np.array([[0.0], [0.0], [0.006]]), t=[0, 1e-7, 2e-7], r=[0]
            )
            restart = runiface(ds, str(directory / 'restart.h5'), quiet=quiet)
            np.testing.assert_allclose(restart.eqsys.f_hot[0], outputs['eigenmode'].eqsys.f_hot[-1], rtol=1e-12)
            np.testing.assert_array_equal(restart.other.hottail.Drr[0], 0)
            np.testing.assert_array_equal(restart.other.hottail.Drr[1], D[-1])
            np.testing.assert_array_equal(restart.other.hottail.Drr_detrapping[:], cap)
            np.testing.assert_array_equal(restart.other.hottail.kappa_detrapping[:], kappa)
            restart.close()
            for output in outputs.values():
                output.close()

            # Archived results stay readable. Running their retired physics is rejected.
            archive = directory / 'eigenmode.h5'
            with h5py.File(archive, 'r+') as raw:
                raw['settings/eqsys/f_hot/transport/detrapping_mode'][...] = 2
            archived = DREAMOutput(str(archive))
            old = archived.settings.eqsys.f_hot.transport
            assert old.detrapping_mode == 2
            try:
                old.verifySettings()
            except Transport.UnsupportedDetrappingSettings:
                pass
            else:
                raise AssertionError('Accepted archived generalized settings.')
            archived.close()

            for mode, scale, message in [(1, 1, 'eigenmode'), (2, 1, 'eigenmode'),
                                          (None, 1, 'eigenmode'), (3, 2, 'detrapping_scale=1')]:
                path = directory / 'legacy_settings.h5'
                settings().save(str(path))
                with h5py.File(path, 'r+') as raw:
                    group = raw['eqsys/f_hot/transport']
                    if mode is None:
                        del group['detrapping_mode']
                    else:
                        group['detrapping_mode'][...] = mode
                    group['detrapping_scale'] = np.array([float(scale)])
                old = DREAMSettings(str(path)).eqsys.f_hot.transport
                try:
                    old.todict()
                except Transport.UnsupportedDetrappingSettings:
                    pass
                else:
                    raise AssertionError(f'Re-serialized retired settings: {mode}, {scale}')
                check_rejected_file(executable, path, message)

            ds = settings()
            ds.hottailgrid.setNxi(3)
            path = directory / 'invalid_pitch.h5'
            ds.save(str(path))
            check_rejected_file(executable, path, 'only applicable in isotropic')

        dreamtests.print_ok('Eigenmode RR settings, coefficients, diagnostics and restart passed.')
        return True
    except Exception as error:
        dreamtests.print_error(f'Detrapping RR regression failed: {error}')
        if not quiet:
            raise
        return False
