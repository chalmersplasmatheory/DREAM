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
is sought. In this case, the ion densities can be controlled (via prescribing
the values or radial transport) by the total ion density :math:`n_i = \sum_j n_i^{(j)}`.

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

When running kinetic simulations, ẗhe density of all ion species, as well as of electrons,
are locally conserved. The net charge of the plasma is also strictly conserved. Therefore,
in order to satisfy quasi-charge neutrality, :math:`n_e = \sum_{ij} Z_{0j}n_i^{(j)}`, the initial plasma must do so.

When setting the initial electron distribution, the ions have been equipped with the function
``getFreeElectronDensity()`` to help the user set the correct density value of the electrons.

Example
^^^^^^^

An example of how to initialize the electron distribution to a Maxwellian with the correct value
for the electron density is provided by the following:

.. code-block:: python

   T_initial = 100 # eV
   n0_initial, rn0 = ds.eqsys.n_i.getFreeElectronDensity()

   ds.eqsys.f_hot.setInitialProfiles(n0=n_initial, rn0=rn0, T0=T_initial)

.. warning::

   The method ``getFreeElectronDensity()`` must be called *after* all initial ion densities have been prescribed.

Ionization models
-----------------

When prescribing initial ion densities as ``DYNAMIC``, their charge states will be evolved
via rate equations including the effects of ionization and recombination. In DREAM, `cold`
electrons will contribute ionization via ``ADAS`` rate coefficients which have been averaged
over a Maxwellian distribution. There is also the option to include a kinetic model of ionization 
due to the fast (non-thermal) electrons, valid up to arbitrarily high energies. The kinetic
ionization model is described in `N A Garland et al, Phys Plasmas 27, 040702 (2020) <https://doi.org/10.1063/5.0003638>`_,
but where DREAM has fitted parameters in the model in order to reproduce the ADAS data as 
accurately as possible.

The ionization model to be used is controlled with the help of three settings:

+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Name                                   | Description                                                                                                                     |
+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONIZATION_MODE_FLUID``              | Only `cold` electrons contribute to ionization, via `ADAS` rate coefficients.                                                   |
+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONIZATION_MODE_KINETIC``            | All electrons in the `f_hot` and `f_re` distributions contribute to ionization, employing the full jacobian.                    |
+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| ``IONIZATION_MODE_KINETIC_APPROX_JAC`` | All electrons in the `f_hot` and `f_re` distributions contribute to ionization, but only the `FLUID` jacobian is used           |
+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------+

In a simulation with many ion species (for example including one hydrogen species and one argon species, 
yielding N=21), the jacobian matrix in mode `KINETIC` will pick up a large number of non-zero elements,
which will generally make simulations extremely heavy. Therefore, it is recommended to include kinetic 
ionization effects via `KINETIC_APPROX_JAC` unless this clearly leads to convergence issues.

Example
^^^^^^^

Kinetic ionization can be activated in a DREAM simulation as follows:

.. code-block:: python

   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC_APPROX_JAC)



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

