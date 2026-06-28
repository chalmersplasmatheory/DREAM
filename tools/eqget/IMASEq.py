# Load equilibrium from IMAS

import imas
import numpy as np
from EQDSK import EQDSK


class IMASEq(EQDSK):
    

    def __init__(
        self, uri, dd_version=None, xml_path=None,
        equilibrium_occurrence=0, wall_occurrence=0, time_index=0,
        profiles_2d_index=0, limiter_unit_index=0,
        override_psilim=False
    ):
        """
        Constructor.

        :param uri: IMAS URI to open in read mode.
        :param dd_version: IMAS data dictionary version.
        :param dd_xml_path: Path to a custom IDSDef.xml.
        """
        dbentry = imas.DBEntry(uri, 'r', dd_version=dd_version, xml_path=xml_path)
        self.eqdsk = None

        try:
            self.eqdsk = self.load_profiles(
                dbentry, equilibrium_occurrence=equilibrium_occurrence,
                wall_occurrence=wall_occurrence, time_index=time_index,
                profiles_2d_index=profiles_2d_index, limiter_unit_index=limiter_unit_index
            )
        finally:
            dbentry.close()

        if self.eqdsk is None:
            raise ValueError("Failed to load equilibrium data from IMAS.")

        super().__init__(self, self.eqdsk, override_psilim=override_psilim)


    def _value(node, *, name):
        if not node.has_value:
            raise ValueError(f"Required IMAS field is missing: {name}")
        return node.value


    def _array(node, *, name):
        return np.asarray(_value(node, name=name), dtype=float)


    def _scalar(node, *, name):
        return float(_value(node, name=name))


    def _string(node, *, default=''):
        if node.has_value:
            return default

        value = node.value
        if isinstance(value, bytes):
            return value.decode()

        return str(value)
    

    def _get_summary_comment(dbentry, occurrence):
        try:
            summar = dbentry.get('summary', occurrence=occurrence)
        except Exception:
            return ''

        return self._string(summary.ids_properties.comment, default='')


    def load_profiles(
        self, dbentry, equilibrium_occurrence=0, wall_occurrence=0,
        time_index=0, profiles_2d_index=0, limiter_unit_index=0
    ):
        """
        Load equilibrium data from the database.
        """
        equilibrium = dbentry.get('equilibrium', occurrence=equilibrium_occurrence)
        wall = dbentry.get('wall', occurrence=wall_occurrence)

        if equilibrium.time_slice.size <= time_index:
            raise IndexError(
                f"Time index {time_index} is out of bounds for equilibrium "+
                f"time slices (size {equilibrium.time_slice.size})."
            )

        ts = equilibrium.time_slice[time_index]
        if ts.profiles_2d.size <= profiles_2d_index:
            raise IndexError(
                f"Profiles 2D index {profiles_2d_index} is out of bounds for "+
                f"time slice {time_index} (size {ts.profiles_2d.size})."
            )

        if wall.description_2d.size == 0:
            raise ValueError("Wall description 2D is empty.")
        if wall.description_2d[0].limiter.unit.size <= limiter_unit_index:
            raise IndexError(
                f"Limiter unit index {limiter_unit_index} is out of bounds "+
                f"for wall description 2D (size "+
                f"{wall.description_2d[0].limiter.unit.size})."
            )

        p1d = ts.profiles_1d
        p2d = ts.profiles_2d[profiles_2d_index]
        boundary = ts.boundary.outline
        global_quantities = ts.global_quantities
        limiter = wall.description_2d[0].limiter.unit[limiter_unit_index].outline

        R = self._array(p2d.grid.dim1, name='equilibrium.time_slice[].profiles_2d[].grid.dim1')
        Z = self._array(p2d.grid.dim2, name='equilibrium.time_slice[].profiles_2d[].grid.dim2')

        psimesh = (
            self._array(p1d.psi_norm, name='equilibrium.time_slice[].profiles_1d.psi_norm')
            if p1d.psi_norm.has_value
            else np.linspace(0.0, 1.0, R.size)
        )

        comment = self._get_summary_comment(dbentry, equilibrium_occurrence)
        if not comment:
            comment = self._string(equilibrium.ids_properties.comment, default='Load from IMAS')

        eqdsk = {
            'stitle': comment,
            'ind1': 1,
            'nx': R.size,
            'ny': Z.size,
            'fpol': self._array(p1d.f, name='equilibrium.time_slice[].profiles_1d.f'),
            'ffprime': self._array(p1d.f_df_dpsi, name='equilibrium.time_slice[].profiles_1d.f_df_dpsi'),
            'p': self._array(p1d.pressure, name='equilibrium.time_slice[].profiles_1d.pressure'),
            'pprime': self._array(p1d.dpressure_dpsi, name='equilibrium.time_slice[].profiles_1d.dpressure_dpsi'),
            'psi': self._array(p2d.psi, name='equilibrium.time_slice[].profiles_2d[].psi'),
            'q': self._array(p1d.q, name='equilibrium.time_slice[].profiles_1d.q'),
            'nbbound': boundary.r.size,
            'nblim': limiter.r.size,
            'rplas': self._array(boundary.r, name='equilibrium.time_slice[].boundary.outline.r'),
            'zplas': self._array(boundary.z, name='equilibrium.time_slice[].boundary.outline.z'),
            'rlim': self._array(limiter.r, name='wall.description_2d[].limiter.unit[].outline.r'),
            'zlim': self._array(limiter.z, name='wall.description_2d[].limiter.unit[].outline.z'),
            'rboxlen': R.max() - R.min(),
            'zboxlen': Z.max() - Z.min(),
            'rcentr': R.min() + 0.5 * (R.max() - R.min()),
            'rleft': R.min(),
            'zmid': Z.min() + 0.5 * (Z.max() - Z.min()),
            'raxis': self._scalar(global_quantities.magnetic_axis.r, name='equilibrium.time_slice[].global_quantities.magnetic_axis.r'),
            'zaxis': self._scalar(global_quantities.magnetic_axis.z, name='equilibrium.time_slice[].global_quantities.magnetic_axis.z'),
            'psiaxis': self._scalar(global_quantities.psi_axis, name='equilibrium.time_slice[].global_quantities.psi_axis'),
            'psiedge': self._scalar(global_quantities.psi_boundary, name='equilibrium.time_slice[].global_quantities.psi_boundary'),
            'bcentr': self._scalar(global_quantities.magnetic_axis.b_field_tor, name='equilibrium.time_slice[].global_quantities.magnetic_axis.b_field_tor'),
            'ip': self._scalar(global_quantities.ip, name='equilibrium.time_slice[].global_quantities.ip'),
            'R': R,
            'Z': Z,
            'psimesh': psimesh,
            'rhopsi': np.sqrt(np.clip(psimesh, 0.0, None)),
            'cocos': 1,
        }

        return eqdsk


