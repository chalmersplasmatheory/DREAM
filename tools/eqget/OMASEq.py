# Load equilibrium from IMAS or OMAS HDF5 via OMAS.

import os
from pathlib import Path
from urllib.parse import unquote

import numpy as np
from omas import load_omas_h5, load_omas_imas

from EQDSK import EQDSK


class OMASEq(EQDSK):

    def __init__(
        self, uri, dd_version=None, backend=None,
        equilibrium_occurrence=0, wall_occurrence=0, time_index=0,
        profiles_2d_index=0, limiter_description_2d_index=0,
        limiter_unit_index=0, consistency_check=False, verbose=False,
        override_psilim=False
    ):
        """
        Constructor.

        :param uri: IMAS URI or path to an OMAS HDF5 file.
        :param dd_version: IMAS data dictionary version to request through OMAS.
        :param backend: Optional backend override. If omitted, the URI backend is used.
        :param override_psilim: Passed to ``EQDSK`` processing.
        """
        ods = self._load_omas(
            uri,
            dd_version=dd_version,
            backend=backend,
            equilibrium_occurrence=equilibrium_occurrence,
            wall_occurrence=wall_occurrence,
            consistency_check=consistency_check,
            verbose=verbose
        )

        eqdsk = self.load_profiles(
            ods,
            time_index=time_index,
            profiles_2d_index=profiles_2d_index,
            limiter_description_2d_index=limiter_description_2d_index,
            limiter_unit_index=limiter_unit_index
        )

        super().__init__(eqdsk, override_psilim=override_psilim)


    @staticmethod
    def _looks_like_h5_file(uri):
        if isinstance(uri, os.PathLike):
            uri = os.fspath(uri)

        if not isinstance(uri, str):
            return False

        path = Path(uri)
        return path.suffix.lower() in ('.h5', '.hdf5')


    @staticmethod
    def _parse_uri(uri, backend=None):
        if not isinstance(uri, str) or not uri.startswith('imas:'):
            raise ValueError(f"Unsupported IMAS URI: {uri!r}")

        scheme, _, query = uri.partition('?')
        parsed_backend = scheme.split(':', 1)[1].strip().upper()
        if not parsed_backend:
            raise ValueError(f"Unable to determine backend from URI: {uri!r}")

        params = {}
        if query:
            for chunk in query.replace('&', ';').split(';'):
                if not chunk:
                    continue

                key, sep, value = chunk.partition('=')
                if not sep:
                    params[unquote(key)] = ''
                else:
                    params[unquote(key)] = unquote(value)

        resolved_backend = (backend or parsed_backend).upper()
        machine = params.get('database', params.get('machine', params.get('path')))
        pulse = params.get('pulse', params.get('shot'))
        run = params.get('run', '0')
        user = params.get('user', os.environ.get('USER', 'dummy_user'))

        if machine is None:
            raise ValueError(
                f"Unable to determine IMAS database/path from URI: {uri!r}"
            )
        if pulse is None:
            raise ValueError(f"Unable to determine IMAS pulse from URI: {uri!r}")

        try:
            pulse = int(pulse)
            run = int(run)
        except ValueError as ex:
            raise ValueError(
                f"Invalid pulse/run in IMAS URI: pulse={pulse!r}, run={run!r}"
            ) from ex

        return {
            'backend': resolved_backend,
            'machine': machine,
            'pulse': pulse,
            'run': run,
            'user': user
        }


    @staticmethod
    def _load_omas(
        uri, dd_version=None, backend=None, equilibrium_occurrence=0,
        wall_occurrence=0, consistency_check=False, verbose=False
    ):
        if OMASEq._looks_like_h5_file(uri):
            return load_omas_h5(
                os.fspath(uri),
                consistency_check=consistency_check,
                imas_version=dd_version
            )

        params = OMASEq._parse_uri(uri, backend=backend)
        occurrence = {
            'equilibrium': equilibrium_occurrence,
            'wall': wall_occurrence
        }

        return load_omas_imas(
            user=params['user'],
            machine=params['machine'],
            pulse=params['pulse'],
            run=params['run'],
            occurrence=occurrence,
            paths=[['equilibrium'], ['wall']],
            imas_version=dd_version,
            consistency_check=consistency_check,
            verbose=verbose,
            backend=params['backend']
        )


    @staticmethod
    def _require(ods, path):
        value = ods.get(path, None)
        if value is None:
            raise ValueError(f"Required OMAS field is missing: {path}")
        return value


    @staticmethod
    def _array(ods, path):
        return np.asarray(OMASEq._require(ods, path), dtype=float)


    @staticmethod
    def _scalar(ods, path):
        return float(OMASEq._require(ods, path))


    @staticmethod
    def _string(ods, path, default=''):
        value = ods.get(path, None)
        if value is None:
            return default
        if isinstance(value, bytes):
            return value.decode()
        return str(value)


    def load_profiles(
        self, ods, time_index=0, profiles_2d_index=0,
        limiter_description_2d_index=0, limiter_unit_index=0
    ):
        """
        Load equilibrium data from an OMAS ODS.
        """
        ts_path = f'equilibrium.time_slice.{time_index}'
        p2d_path = f'{ts_path}.profiles_2d.{profiles_2d_index}'
        limiter_path = (
            'wall.description_2d.'
            f'{limiter_description_2d_index}.limiter.unit.{limiter_unit_index}.outline'
        )

        R = self._array(ods, f'{p2d_path}.grid.dim1')
        Z = self._array(ods, f'{p2d_path}.grid.dim2')
        rplas = self._array(ods, f'{ts_path}.boundary.outline.r')
        zplas = self._array(ods, f'{ts_path}.boundary.outline.z')
        if ods.get(f'{limiter_path}.r', None) is None or ods.get(f'{limiter_path}.z', None) is None:
            rlim = rplas.copy()
            zlim = zplas.copy()
        else:
            rlim = self._array(ods, f'{limiter_path}.r')
            zlim = self._array(ods, f'{limiter_path}.z')

        psimesh_path = f'{ts_path}.profiles_1d.psi_norm'
        if ods.get(psimesh_path, None) is None:
            psimesh = np.linspace(0.0, 1.0, R.size)
        else:
            psimesh = np.asarray(ods[psimesh_path], dtype=float)

        comment = self._string(ods, 'summary.ids_properties.comment')
        if not comment:
            comment = self._string(
                ods, 'equilibrium.ids_properties.comment', default='Load from OMAS'
            )

        return {
            'stitle': comment,
            'ind1': 1,
            'nx': R.size,
            'ny': Z.size,
            'fpol': self._array(ods, f'{ts_path}.profiles_1d.f'),
            'ffprime': self._array(ods, f'{ts_path}.profiles_1d.f_df_dpsi'),
            'p': self._array(ods, f'{ts_path}.profiles_1d.pressure'),
            'pprime': self._array(ods, f'{ts_path}.profiles_1d.dpressure_dpsi'),
            'psi': self._array(ods, f'{p2d_path}.psi'),
            'q': self._array(ods, f'{ts_path}.profiles_1d.q'),
            'nbbound': rplas.size,
            'nblim': rlim.size,
            'rplas': rplas,
            'zplas': zplas,
            'rlim': rlim,
            'zlim': zlim,
            'rboxlen': R.max() - R.min(),
            'zboxlen': Z.max() - Z.min(),
            'rcentr': R.min() + 0.5 * (R.max() - R.min()),
            'rleft': R.min(),
            'zmid': Z.min() + 0.5 * (Z.max() - Z.min()),
            'raxis': self._scalar(
                ods, f'{ts_path}.global_quantities.magnetic_axis.r'
            ),
            'zaxis': self._scalar(
                ods, f'{ts_path}.global_quantities.magnetic_axis.z'
            ),
            'psiaxis': self._scalar(ods, f'{ts_path}.global_quantities.psi_axis'),
            'psiedge': self._scalar(ods, f'{ts_path}.global_quantities.psi_boundary'),
            'bcentr': self._scalar(
                ods, f'{ts_path}.global_quantities.magnetic_axis.b_field_tor'
            ),
            'ip': self._scalar(ods, f'{ts_path}.global_quantities.ip'),
            'R': R,
            'Z': Z,
            'psimesh': psimesh,
            'rhopsi': np.sqrt(np.clip(psimesh, 0.0, None)),
            'cocos': 1,
        }
