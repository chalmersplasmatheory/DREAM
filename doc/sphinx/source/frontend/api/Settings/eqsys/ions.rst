.. _ds-eqsys-ions:

Ions
====
Ions play an important role in runaway electron simulations. They are often the
main source of pitch angle scattering, cause electrons to lose energy through
the production of bremsstrahlung and provide a source of free electrons as they
ionize. The ions also affect the energy balance of the plasma as they further
ionize and recombine, which is particularly important during tokamak
disruptions.

In DREAM, the user can specify an arbitrary number of ions to include in the
simulation. During the simulation DREAM also keeps track of the charge states
of the ions which enters into calculations of the plasma energy balance, as well
as the interaction of electrons with partially ionized atoms.

.. warning::

   Since ions and electrons are specified separately (the latter by prescribing
   either the electron distribution function or cold electron density), it is
   possible to set up a simulation in which ions and electrons are unbalanced.
   In kinetic mode, DREAM will always preserve quasi-neutrality by adding as
   many cold electrons as necessary, which may not be what the user desires.

   To ensure that the number of electrons and ions agree, use the
   ``getFreeElectronDensity()`` method in the ``Ions`` class after specifying
   the ions for the simulation, and use this density to initialize the electrons
   in the simulation.

.. contents:: Page overview
   :local:
   :depth: 3

Ion models
----------
DREAM provides three different modes for evolving ions. It is possible to evolve
different ions species with different modes in the same simulation:

(1) Prescribe the evolution in time and space of the ion densities and charge states.
(2) Assume that ion charge states are distributed as to be in equilibrium in every time step of the simulation.
(3) Evolve the ion charge states using an ion rate equation.

In mode (1) it is up to the user to provide the full evolution of the ions.
Mode (3) evolves the ion charge states by solving the ion rate equation

.. math::

   \frac{\partial n_i^{(j)}}{\partial t} =
       \left( I_i^{(j-1)} n_{\rm cold} + \mathcal{I}_i^{(j-1)} \right) n_i^{(j-1)} -
       \left( I_i^{(j)} n_{\rm cold} + \mathcal{I}_i^{(j)} \right) n_i^{(j)}\\
       + R_i^{(j+1)} n_i^{(j+1)} n_{\rm cold} - R_i^{(j)} n_i^{(j)} n_{\rm cold}

where :math:`n_i^{(j)}` denotes the density of ion species :math:`i` in charge
state :math:`j`, :math:`I_i^{(j)}` is the rate at which ion species :math:`i`
ionizes from charge state :math:`j` to charge state :math:`j+1`,
:math:`R_i^{(j)}` is the rate at which ion species :math:`i` recombines from
charge state :math:`j` to charge state :math:`j-1`, and
:math:`\mathcal{I}_i^{(j)}` is the rate at which ion species :math:`i` in charge
state :math:`j` is ionized due to collisions with fast electrons. The ionization
and recombination rates :math:`I_i^{(j)}` and :math:`R_i^{(j)}` are taken from
the `OPEN-ADAS database <https://open.adas.ac.uk/>`.

In ion mode (2), the above ion rate equation is also solved, but with
:math:`\partial n_i^{(j)} / \partial t = 0` such that the equilibrium solution
is sought.

The available ion modes are:

+----------------------+-----------------------------------------------------------------------+
| Name                 | Description                                                           |
+======================+=======================================================================+
| ``IONS_PRESCRIBED``  | The ion densities and charge states are prescribed in time and space. |
+----------------------+-----------------------------------------------------------------------+
| ``IONS_EQUILIBRIUM`` | Ion charge states are assumed to be in equilibrium.                   |
+----------------------+-----------------------------------------------------------------------+
| ``IONS_DYNAMIC``     | Ion charge states are evolved according to an ion rate equation.      |
+----------------------+-----------------------------------------------------------------------+

The above modes are the only three models available for evolving ions in DREAM.
However, to simplify initialization of ion densities, the following pseudo-modes
are also available when adding ions to the simulation:

+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Name                               | Description                                                                                                                     |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONS_DYNAMIC_NEUTRAL``           | Ion densities are evolved dynamically, with all ions located in the neutral charge state initially.                             |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONS_DYNAMIC_FULLY_IONIZED``     | Ion densities are evolved dynamically, with all ions located in the fully ionized charge state initially.                       |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONS_PRESCRIBED_NEUTRAL``        | Ion densities are prescribed, with all ions located in the neutral charge state.                                                |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONS_PRESCRIBED_FULL_IONIZED``   | Ion densities are prescribed, with all ions located in the fully ionized charge state.                                          |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONS_EQUILIBRIUM_NEUTRAL``       | Ion densities are assumed in equilibrium at every time step, with all ions located in the neutral charge state initially.       |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONS_EQUILIBRIUM_FULLY_IONIZED`` | Ion densities are evolved in equilibrium at every time step, with all ions located in the fully ionized charge state initially. |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+

ADAS and NIST
-------------

.. todo::

   Describe how to use the ``tools/get_adas.py`` and ``tools/get_nist.py``
   scripts.

Adding ions
-----------

.. todo::

   Describe how to add ions to simulations.

Free electron density
---------------------

.. todo::

   Describe how and when to use the ``getFreeElectronDensity()`` method.

Ionization models
-----------------

.. todo::

   Describe the different ionization models available.

Tritium
-------

.. todo::

   Describe the specification of tritium.

Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.Ions.Ions
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

.. autoclass:: DREAM.Settings.Equations.IonSpecies.IonSpecies
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

