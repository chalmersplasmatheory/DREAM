.. _ds-eqsys-molecular-reactions:

Molecular reactions
===================
This page lists all molecular reactions which are available in DREAM. The
reactions are defined in the Python module
``py/DREAM/Settings/MolecularReactionData.py``.

For instructions on how to add molecules to a simulation and how to enable
reactions, see :ref:`Molecules<ds-eqsys-ions-molecules>` on the
:ref:`ions<ds-eqsys-ions>` page.

.. contents:: Page overview
   :local:
   :depth: 2

Reading the tables
------------------
Each reaction is listed with

- **Name**: the unique name of the reaction. This is the name used when
  replacing a built-in rate with ``setReactionRate()``, and it is listed in
  ``ds.eqsys.n_i.enabled_reaction_names`` once the reaction has been enabled.
- **Reaction**: the reactants and products. The electron :math:`\mathrm{e}^-`
  refers to the cold electron population (``n_cold``), while
  :math:`\mathrm{e}^-_\mathrm{RE}` refers to the runaway electrons (``n_re``).
- **Classification**: the group to which the reaction belongs. This is the key
  used when selecting reactions with ``enableReactions()``.

In the reaction rate equation, each reaction :math:`r` with reactants
:math:`A` and :math:`B` contributes a reaction rate per unit volume

.. math::

   \mathcal{R}_r = \langle\sigma v\rangle_r n_A n_B,

where :math:`\langle\sigma v\rangle_r` is the rate coefficient (in
:math:`\mathrm{m}^3/\mathrm{s}`) and :math:`n_A`, :math:`n_B` are the reactant
densities (in :math:`\mathrm{m}^{-3}`), so that :math:`\mathcal{R}_r` is given
in :math:`\mathrm{m}^{-3}\mathrm{s}^{-1}`. Each reaction removes the reactants
and adds the products at this rate in the ion rate equation.

Ionization
----------
Electron-impact ionization of molecules. Molecules are not included in the ADAS
database, so these reactions replace the ADAS ionization used for atomic
species.

.. molecularreactionlist::
   :classification: IONIZATION

Recombination
-------------
Dissociative recombination of molecular ions with cold electrons.

.. molecularreactionlist::
   :classification: RECOMBINATION

Dissociation
------------
Electron-impact dissociation of molecules and molecular ions. Several reactions
can have the same reactants but different products (branches). Selecting the
reactants in ``enableReactions()`` enables all branches.

.. molecularreactionlist::
   :classification: DISSOCIATION

Heavy-particle dissociation
---------------------------
Reactions between two heavy particles where molecules are formed or broken up
(``DISSOCIATION_MULTIPLE``).

.. molecularreactionlist::
   :classification: DISSOCIATION_MULTIPLE

Charge exchange
---------------
Transfer of an electron between two different species, or between different
charge states of the same species.

.. molecularreactionlist::
   :classification: CHARGE_EXCHANGE

Resonant charge exchange
------------------------
Charge exchange between the neutral and the ion of the same species. These
reactions do not change any densities, and are therefore not included in the
ion rate equation.

.. molecularreactionlist::
   :classification: CHARGE_EXCHANGE_RESONANT

Runaway ionization
------------------
Ionization of molecules by runaway electrons. These reactions are only
included when the ionization mode is ``IONIZATION_MODE_FLUID_RE``, and are
rescaled by the same factor as the atomic runaway ionization term. For atomic
species the corresponding term is ``IonFluidRunawayIonizationTerm``.

.. molecularreactionlist::
   :classification: RUNAWAY

Adding a new reaction
---------------------
New reactions are added as entries in the list ``MOLECULAR_REACTIONS`` in
``py/DREAM/Settings/MolecularReactionData.py``, with any numerical
coefficients placed in ``MolecularReactionDataCoefficients.py``. The rates are
converted into C++ tables by ``tools/MolecularRates/generate_molecular_rates.py``
when DREAM is built, so DREAM must be recompiled after a reaction is added.
The tables on this page are updated automatically the next time the
documentation is built.

A reaction entry has the following form:

.. code-block:: python

   {
       "name": "D_2_ionization",
       "reactants": [("e", -1, 1), ("D2", 0, 1)],          # (name, Z0, number of particles)
       "products": [("e", -1, 2), ("D2", 1, 1)],
       "classification": "IONIZATION",
       "data_type": "AMJUEL_POLYNOMIAL",                     # How the rate is evaluated
       "data": {...},                                       # Coefficients for 'data_type'
       "temperature_input": {"kind": "electron"},            # Temperature used in the rate [eV]
       "density_input": {"kind": "electron"}                 # Density used in the rate [m^-3]
   }

.. note::

   If a reaction involves a species which is not yet used by any reaction, its
   charge ``Z`` must be added to ``SPECIES_Z`` in
   ``MolecularReactionData.py``. If the species is a molecule, it must also be
   added to the table in ``src/MoleculeHandler.cpp``, which gives its mass and
   the charge of its heaviest nucleus.
