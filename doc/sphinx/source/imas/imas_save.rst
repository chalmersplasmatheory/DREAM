.. _imas_save:

Exporting DREAM output to IMAS
==============================

DREAM output files can be converted into IMAS format using the scripts located
under ``tools/h5_to_IMAS/``. The main conversion script is
``h5file_to_imas.py``, which reads one DREAM output file and writes a
corresponding IMAS entry.

.. contents:: Page overview
   :local:
   :depth: 3


Supported IDSs
--------------

Both exporters support the same set of IDSs, selected with ``--ids``:

* ``plasma_profiles``
* ``runaway_electrons``
* ``spi``
* ``equilibrium``
* ``radiation``
* ``summary``

The scripts also generate a Markdown mapping report. This report contains the
exact node-by-node mapping for the exported file or simulation, including any
skipped nodes, missing sources and warnings. The example reported in
`PR #524 <https://github.com/chalmersplasmatheory/DREAM/pull/524>`_ writes
140 unique mappings across the six IDSs above.

.. list-table:: Implemented DREAM to IMAS mapping overview
   :header-rows: 1
   :widths: 18 42 40

   * - IDS
     - Main IMAS content written
     - Main DREAM sources
   * - ``plasma_profiles``
     - Top-level time and vacuum toroidal field; 1-D electron profiles;
       ``j_total`` and ``j_ohmic``; conductivity; ``zeff``; radial grids; ion
       and neutral composition; ion and neutral charge-state densities.
     - ``/grid/t``, ``/grid/B0``, ``/grid/R0``, ``/eqsys/T_cold``,
       ``/eqsys/n_tot`` or ``/eqsys/n_cold``, ``/eqsys/n_i``,
       ``/eqsys/j_tot``, ``/eqsys/j_ohm``, ``/eqsys/E_field``,
       ``/other/fluid/Zeff``, ``/other/fluid/conductivity`` and grid geometry
       datasets.
   * - ``runaway_electrons``
     - Top-level time and vacuum toroidal field; runaway density and current
       density; Dreicer, hot-tail, Compton, tritium and total source terms;
       critical electric fields and critical momenta; 1-D radial grids; total
       runaway current.
     - ``/eqsys/n_re``, ``/eqsys/j_re``, ``/other/fluid/runawayRate``,
       ``/other/fluid/gammaDreicer``, ``/other/fluid/gammaHottail``,
       ``/other/fluid/gammaCompton``, ``/other/fluid/gammaTritium``,
       ``/other/fluid/EDreic``, ``/other/fluid/Ectot``,
       ``/other/fluid/pCrit``, ``/other/fluid/pCritHottail`` and grid
       geometry datasets.
   * - ``spi``
     - Injector metadata; pellet composition; shard positions, velocities and
       volumes; fragment mass-centre velocities.
     - ``/settings/eqsys/spi/init/*``, ``/settings/eqsys/spi/ZsDrift``,
       ``/settings/eqsys/spi/isotopesDrift``,
       ``/settings/eqsys/n_i/SPIMolarFraction``, ``/eqsys/x_p``,
       ``/eqsys/v_p`` and ``/eqsys/Y_p``.
   * - ``equilibrium``
     - Top-level time and vacuum toroidal field; boundary outline; geometric
       axis; plasma current; poloidal flux references; 1-D equilibrium
       profiles; 2-D poloidal-plane grids with ``r``, ``z`` and ``psi``.
     - ``/grid/t``, ``/grid/B0``, ``/grid/R0``, ``/eqsys/I_p``,
       ``/eqsys/psi_p``, ``/eqsys/j_tot``, ``/grid/eq/*`` and
       ``/grid/geometry/*``.
   * - ``radiation``
     - A single electron-radiation process with global radiated power and 1-D
       emissivity and cumulative power profiles.
     - ``/other/fluid/Tcold_radiation`` aligned to ``/grid/t``, together with
       ``/grid/dr``, ``/grid/VpVol``, ``/grid/R0`` and toroidal-flux geometry
       data.
   * - ``summary``
     - Top-level time; code provenance; global current and thermal-energy
       traces; magnetic-axis and separatrix values; volume averages; runaway
       current and particle content.
     - ``/grid/t``, ``/code/*``, ``/eqsys/I_p``, ``/eqsys/j_ohm``,
       ``/eqsys/j_re``, ``/eqsys/W_cold``, ``/eqsys/W_i``,
       ``/eqsys/n_tot`` or ``/eqsys/n_cold``, ``/eqsys/n_re``,
       ``/eqsys/T_cold``, ``/eqsys/E_field``, ``/other/fluid/Zeff`` and grid
       volume factors.


Single-file export
------------------

Use ``tools/h5_to_IMAS/h5file_to_imas.py`` to convert one DREAM output file:

.. code-block:: bash

   python tools/h5_to_IMAS/h5file_to_imas.py output_CQ_S6.h5

By default this writes to ``dream_imas.nc`` and generates a mapping report next
to the output. Typical options are:

* ``--uri``: output URI passed to ``imas.DBEntry``.
* ``--ids``: select one or more of the supported IDSs.
* ``--dd-version``: request a specific IMAS Data Dictionary version.
* ``--report``: set the mapping-report filename.
* ``--dry-run``: build the IDS objects and mapping report without calling
  ``DBEntry.put()``.

Examples:

.. code-block:: bash

   python tools/h5_to_IMAS/h5file_to_imas.py output_CQ_S6.h5 \
       --uri 'imas:hdf5?path=./dream_imas'

   python tools/h5_to_IMAS/h5file_to_imas.py output_CQ_S6.h5 \
       --ids plasma_profiles runaway_electrons summary \
       --report output_CQ_S6.mapping_report.md

   python tools/h5_to_IMAS/h5file_to_imas.py output_CQ_S6.h5 --dry-run

The mapping report is useful even for dry runs, because it shows which IMAS
nodes were filled, which derivations were used and which nodes were skipped by
the selected Data Dictionary.


Multi-stage simulation export
-----------------------------

Use ``tools/h5_to_IMAS/export_multiple_sims_to_IMAS.py`` when each simulation
is split across several DREAM output files:

.. code-block:: bash

   python tools/h5_to_IMAS/export_multiple_sims_to_IMAS.py /path/to/simulations \
       --pulse-start 120000 \
       --pulse-end 129999 \
       --write-simulation-reports

The exporter looks for files matching ``output*.h5`` and treats each simulation
directory as one IMAS entry. It first searches for stage files in ``output/``
subdirectories, then falls back to files directly under the root path or one
level below it.

Stage handling
^^^^^^^^^^^^^^

The stage files are ordered by filename pattern:

* ``init``
* ``injection_start``
* ``injection_continuation``
* ``injection``
* ``cq``
* everything else, sorted naturally and warned about

Files classified as ``init`` are ignored by the batch exporter. For every later
stage, the script:

1. reads the local DREAM time grid from ``/grid/t``,
2. offsets the stage so that it starts after the previous stage,
3. trims overlapping samples,
4. writes the first slice with ``DBEntry.put()``,
5. appends later slices with ``DBEntry.put_slice()``.

The ``spi`` IDS is handled slightly differently: the script concatenates the
fragment data in memory and writes one combined IDS with ``DBEntry.put()``.

The most useful options are:

* ``--stage-glob``: glob used to find the stage files.
* ``--uri-template``: URI template for each exported simulation. The default is
  ``imas:hdf5?user={user};pulse={pulse};run={run};database={database};version={version}``.
* ``--imas-user``, ``--database``, ``--version``: fill the default URI
  template.
* ``--pulse-start``, ``--pulse-end`` and ``--run``: pulse/run numbers assigned
  to exported simulations.
* ``--ids`` and ``--dd-version``: same meaning as for the single-file export.
* ``--report``: batch report path.
* ``--uri-map``: TSV file listing the output URI for each simulation folder.
* ``--write-simulation-reports``: also write one mapping report in each
  simulation directory.
* ``--dry-run``: convert the stages and write reports, but do not create IMAS
  entries.
* ``--profile``: print and record coarse timing information.


Full mapping reference
----------------------

The table below lists all the mappings from DREAM to IMAS.

+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| **IMAS node**                                                         | **Units**    | **DREAM source / derivation**                                                                                    |
+=======================================================================+==============+==================================================================================================================+
| ``equilibrium.time_slice.global_quantities.psi_magnetic_axis``        | Wb           | ``/eqsys/psi_p[:,*] * 2pi``                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.boundary.outline.r``                         | m            | ``/grid/eq/RMinusR0_f[:,-1] + R0``                                                                               |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.boundary.outline.z``                         | m            | ``/grid/eq/ZMinusZ0_f[:,-1] + Z0``                                                                               |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.boundary.geometric_axis.r``                  | m            | derived from ``/grid/eq/RMinusR0_f``                                                                             |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.boundary.geometric_axis.z``                  | m            | derived from ``/grid/eq/ZMinusZ0_f``                                                                             |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.global_quantities.ip``                       | A            | ``/eqsys/I_p``                                                                                                   |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.global_quantities.psi_boundary``             | Wb           | ``/eqsys/psi_p[:,-1] * 2pi``                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.global_quantities.psi_magnetic_axis``        | Wb           | ``/eqsys/psi_p[:,*] * 2pi``                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.b_field_average``                | T            | derived from ``/grid/geometry/Bmin`` and ``/grid/geometry/FSA_BOverBmin``                                        |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.b_field_max``                    | T            | ``/grid/geometry/Bmax``                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.b_field_min``                    | T            | ``/grid/geometry/Bmin``                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.dvolume_dpsi``                   | m^3.Wb^-1    | derived                                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.gm1``                            | m^-2         | derived from ``/grid/geometry/FSA_R02OverR2`` and ``/grid/R0``                                                   |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.gm5``                            | T^2          | derived from ``/grid/geometry/Bmin`` and ``/grid/geometry/FSA_BOverBmin2``                                       |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.j_phi``                          | A.m^-2       | ``/eqsys/j_tot``                                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.phi``                            | Wb           | derived from ``phi_tor``                                                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.psi``                            | Wb           | ``/eqsys/psi_p * 2pi``                                                                                           |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.psi_norm``                       | 1            | derived                                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.r_outboard``                     | m            | derived from R0 and r                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.rho_tor``                        | m            | derived from ``/grid/geometry/toroidalFlux``                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.rho_tor_norm``                   | 1            | derived from ``rho_tor``                                                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_1d.trapped_fraction``               | 1            | derived from ``/grid/geometry/effectivePassingFraction``                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_2d.grid.dim1``                      | mixed        | ``/eqsys/psi_p * 2pi``                                                                                           |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_2d.grid.dim2``                      | mixed        | ``/grid/eq/theta``                                                                                               |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_2d.grid_type.description``          | n/a          | ``poloidal_plane_coordinates_identifier``                                                                        |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_2d.psi``                            | Wb           | ``/grid/eq/RMinusR0_f``, ``/grid/eq/ZMinusZ0_f``; cell-edge psi derived from ``/eqsys/psi_p``                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_2d.r``                              | m            | ``/grid/eq/RMinusR0_f``, ``/grid/eq/ZMinusZ0_f``; cell-edge psi derived from ``/eqsys/psi_p``                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.profiles_2d.z``                              | m            | ``/grid/eq/RMinusR0_f``, ``/grid/eq/ZMinusZ0_f``; cell-edge psi derived from ``/eqsys/psi_p``                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time_slice.time``                                       | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.time``                                                  | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.vacuum_toroidal_field.b0``                              | T            | ``/grid/B0``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``equilibrium.vacuum_toroidal_field.r0``                              | m            | ``/grid/R0``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.ion.density``                           | m^-3         | ``/eqsys/n_i`` summed over charged states :math:`Z_0=1\ldots Z`                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.ion.element.a``                         | u            | ``/settings/eqsys/n_i/Z``, ``/settings/eqsys/n_i/isotopes``                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.ion.element.z_n``                       | e            | ``/settings/eqsys/n_i/Z``, ``/settings/eqsys/n_i/isotopes``                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.neutral.density``                       | m^-3         | ``/eqsys/n_i`` for :math:`Z_0=0`                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.neutral.element.a``                     | u            | ``/settings/eqsys/n_i/Z``, ``/settings/eqsys/n_i/isotopes``                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.neutral.element.z_n``                   | e            | ``/settings/eqsys/n_i/Z``, ``/settings/eqsys/n_i/isotopes``                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.neutral.state.density``                 | m^-3         | ``/eqsys/n_i`` for :math:`Z0=0`                                                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.conductivity_parallel``                 | ohm^-1.m^-1  | ``/other/fluid/conductivity``                                                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.e_field.parallel``                      | V.m^-1       | ``/eqsys/E_field``                                                                                               |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.electrons.density``                     | m^-3         | ``/eqsys/n_tot``                                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.electrons.density_thermal``             | m^-3         | ``/eqsys/n_cold``                                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.electrons.temperature``                 | eV           | ``/eqsys/T_cold``                                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.grid.psi``                              | Wb           | ``/eqsys/psi_p``                                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.grid.psi_boundary``                     | Wb           | ``/eqsys/psi_p[:,-1] * 2pi``                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.grid.psi_magnetic_axis``                | Wb           | ``/eqsys/psi_p[:,*] * 2pi``                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.grid.rho_pol_norm``                     | 1            | derived                                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.grid.rho_tor``                          | m            | ``/grid/geometry/toroidalFlux``                                                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.grid.rho_tor_norm``                     | 1            | derived from ``rho_tor``                                                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.j_ohmic``                               | A.m^-2       | ``/eqsys/j_ohm``                                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.j_total``                               | A.m^-2       | ``/eqsys/j_tot``                                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.time``                                  | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.profiles_1d.zeff``                                  | 1            | ``/other/fluid/Zeff``                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.time``                                              | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.vacuum_toroidal_field.b0``                          | T            | ``/grid/B0``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``plasma_profiles.vacuum_toroidal_field.r0``                          | m            | ``/grid/R0``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.global_quantities.inside_vessel.power_electrons`` |  W           |  volume integral of ``/other/fluid/Tcold_radiation`` using ``dV = /grid/dr * /grid/VpVol * /grid/R0``            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.global_quantities.time``                          | s            | ``/grid/t`` aligned to ``/other/fluid/Tcold_radiation``                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.profiles_1d.electrons.emissivity``                | W.m^-3       | ``/other/fluid/Tcold_radiation``                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.profiles_1d.electrons.power_inside``              | W            | cumulative volume integral of ``/other/fluid/Tcold_radiation`` using ``dV = /grid/dr * /grid/VpVol * /grid/R0``  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.profiles_1d.grid.rho_tor``                        | m            | ``rho_tor = sqrt(/grid/geometry/toroidalFlux/(pi*B0))``                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.profiles_1d.grid.rho_tor_norm``                   | 1            | ``rho_tor_norm = normalized sqrt(/grid/geometry/toroidalFlux/(pi*B0))``                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.profiles_1d.time``                                | s            | ``/grid/t`` aligned to ``/other/fluid/Tcold_radiation``                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.identifier.description``                          | n/a          | ``/other/fluid/Tcold_radiation``                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.identifier.index``                                | 1            | private process identifier for DREAM total electron radiation                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.process.identifier.name``                                 | n/a          | private process identifier for DREAM total electron radiation                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``radiation.time``                                                    | s            | ``/grid/t`` aligned to ``/other/fluid/Tcold_radiation``                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.current_density``                     | A.m^-2       | ``/eqsys/j_re``                                                                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.ddensity_dt_compton``                 | m^-3.s^-1    | ``/other/fluid/gammaCompton``                                                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.ddensity_dt_dreicer``                 | m^-3.s^-1    | ``/other/fluid/gammaDreicer``                                                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.ddensity_dt_hot_tail``                | m^-3.s^-1    | ``/other/fluid/gammaHottail``                                                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.ddensity_dt_total``                   | m^-3.s^-1    | ``/other/fluid/runawayRate``                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.ddensity_dt_tritium``                 | m^-3.s^-1    | ``/other/fluid/gammaTritium``                                                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.density``                             | m^-3         | ``/eqsys/n_re``                                                                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.e_field_critical``                    | V.m^-1       | ``/other/fluid/Ectot``                                                                                           |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.e_field_dreicer``                     | V.m^-1       | ``/other/fluid/EDreic``                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.grid.psi``                            | Wb           | ``/eqsys/psi_p*2pi``                                                                                             |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.grid.psi_boundary``                   | Wb           | ``/eqsys/psi_p[:,-1] * 2pi``                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.grid.psi_magnetic_axis``              | Wb           | ``/eqsys/psi_p[:,*] * 2pi``                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.grid.rho_pol_norm``                   | 1            | derived                                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.grid.rho_tor``                        | m            | ``/grid/geometry/toroidalFlux``                                                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.grid.rho_tor_norm``                   | 1            | derived from ``rho_tor``                                                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.momentum_critical_avalanche``         | kg.m^-1.s^-1 | ``/other/fluid/pCrit * m_e*c``                                                                                   |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.momentum_critical_hot_tail``          | kg.m^-1.s^-1 | ``/other/fluid/pCritHottail * m_e*c``                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.profiles_1d.time``                                |  s           |  ``/grid/t``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.global_quantities.current_phi``                   | A            | derived from ``j_re``                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.time``                                            | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.vacuum_toroidal_field.b0``                        | T            | ``/grid/B0``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``runaway_electrons.vacuum_toroidal_field.r0``                        | m            | ``/grid/R0``                                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.fragment.position.r``                                  | m            | ``/eqsys/x_p[:,*,:] + R0``                                                                                       |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.fragment.position.z``                                  | m            | ``/eqsys/x_p[:,*,:]``                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.fragment.velocity_r``                                  | m.s^-1       | ``/eqsys/v_p[:,*,:]``                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.fragment.velocity_z``                                  | m.s^-1       | ``/eqsys/v_p[:,*,:]``                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.fragment.volume``                                      | m^3          | ``/eqsys/Y_p[:,*]``                                                                                              |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.pellet.core.species.a``                                | u            | ``/settings/eqsys/n_i/names``, ``/settings/eqsys/spi/init/Ninj``                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.pellet.core.species.density``                          | m^-3         | ``/settings/eqsys/n_i/SPIMolarFraction``, ``/eqsys/Y_p`` or ``/settings/eqsys/spi/init/rp``                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.pellet.core.species.name``                             | n/a          | ``/settings/eqsys/n_i/names``, ``/settings/eqsys/spi/init/Ninj``                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.pellet.core.species.z_n``                              | e            | ``/settings/eqsys/n_i/names``, ``/settings/eqsys/spi/init/Ninj``                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.description``                                          | n/a          | DREAM SPI                                                                                                        |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.name``                                                 | n/a          | DREAM SPI                                                                                                        |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.pellet.core.atoms_n``                                  | 1            | ``/settings/eqsys/spi/init/Ninj``                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.velocity_mass_centre_fragments_r``                     | m.s^-1       | ``/eqsys/v_p[*]``                                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.injector.velocity_mass_centre_fragments_z``                     | m.s^-1       | ``/eqsys/v_p[*]``                                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``spi.time``                                                          | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.code.commit``                                               | n/a          | ``/code/commit``                                                                                                 |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.code.description``                                          | n/a          | DREAM code description                                                                                           |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.code.name``                                                 | n/a          | DREAM code name                                                                                                  |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.code.parameters``                                           | n/a          | ``/code/changes``, ``/code/commit``, ``/code/datetime_commit``, ``/code/datetime_simulation``, ``/code/refspec`` |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.code.repository``                                           | n/a          | DREAM public repository                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.code.version``                                              | n/a          | ``/code/refspec``                                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.global_quantities.current_ohm.value``                       | A            | derived from ``/eqsys/j_ohm``                                                                                    |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.global_quantities.energy_electrons_thermal.value``          | J            | volume integral of ``/eqsys/W_cold``                                                                             |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.global_quantities.energy_ion_total_thermal.value``          | J            | volume integral of ``/eqsys/W_i``                                                                                |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.global_quantities.energy_thermal.value``                    | J            | volume integral of ``/eqsys/W_cold`` plus ``/eqsys/W_i`` when available                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.global_quantities.ip.value``                                | A            | ``/eqsys/I_p``                                                                                                   |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.magnetic_axis.e_field_parallel.value``                | V.m^-1       | ``/eqsys/E_field[:,*]``                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.magnetic_axis.n_e.value``                             | m^-3         | ``/eqsys/n_tot[:,*]``                                                                                            |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.magnetic_axis.t_e.value``                             | eV           | ``/eqsys/T_cold[:,*]``                                                                                           |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.magnetic_axis.zeff.value``                            | 1            | ``/other/fluid/Zeff[:,*]``                                                                                       |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.separatrix.e_field_parallel.value``                   | V.m^-1       | ``/eqsys/E_field[:,-1]``                                                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.separatrix.n_e.value``                                | m^-3         | ``/eqsys/n_tot[:,-1]``                                                                                           |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.separatrix.t_e.value``                                | eV           | ``/eqsys/T_cold[:,-1]``                                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.local.separatrix.zeff.value``                               | 1            | ``/other/fluid/Zeff[:,-1]``                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.runaways.current.value``                                    | A            | derived from ``/eqsys/j_re``                                                                                     |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.runaways.current_phi_max.value``                            | A            | maximum absolute derived runaway current                                                                         |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.runaways.particles.value``                                  | 1            | volume integral of ``/eqsys/n_re``                                                                               |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.time``                                                      | s            | ``/grid/t``                                                                                                      |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.volume_average.n_e.value``                                  | m^-3         | volume average of ``/eqsys/n_tot``                                                                               |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.volume_average.t_e.value``                                  | eV           | volume average of ``/eqsys/T_cold``                                                                              |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+
| ``summary.volume_average.zeff.value``                                 | 1            | volume average of ``/other/fluid/Zeff``                                                                          |
+-----------------------------------------------------------------------+--------------+------------------------------------------------------------------------------------------------------------------+

