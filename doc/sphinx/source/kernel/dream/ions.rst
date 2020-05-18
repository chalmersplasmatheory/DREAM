Ions
====
DREAM supports the inclusion of an arbitrary number of ions species. Since we
solve for the density of each ion charge state, each ion contributes
:math:`Z+1` independent fluid (i.e. radially varying) variables to the equation
system. This suggests that a special representation should be used for ions in
the code, in order to keep similar density variables close in memory. As a
result, DREAM combines *all* ion species into a single ion density, ``n_i``.
To keep track of and work with this ion density, a number of helper classes
are introduced, and this section of the documentation describes those.

While the details of the implementation are discussed below, it can be worth
just stating what actually constitutes an ion in DREAM. To add a new ion to the
simulation, the user must provide three things:

- The **name** of the ion species (only used for communicating with the user)
- The **atomic charge** :math:`Z` of the ion species
- The radial **density profile** of every charge state of the ion species

In addition to these, the user will also specify how to evolve the ions during
the simulation. DREAM supports three methods for evolving ions in time, namely

- Prescribed time evolution
- Coronal equilibrium
- Dynamic ionization/recombination balance

The first of these requires the user to provide the full spatiotemporal
evolution of every charge state of the ion species to the code.

The second and third option represent variations of solving the ion rate
equation

.. math::

   \frac{\partial n_i^{(j)}}{\partial t} &=
       \left( I_i^{(j-1)} n_{\rm cold} + \mathcal{I}_i^{(j-1)} \right) n_i^{(j-1)} -
       \left( I_i^{(j)} n_{\rm cold}+ \mathcal{I}_i^{(j)} \right) n_i^{(j)} +\\
       &+ R_i^{(j+1)} n_i^{(j+1)} n_{\rm cold} - R_i^{(j)} n_i^{(j)} n_{\rm cold},

where :math:`n_i^{(j)}` denotes the radial density of charge state :math:`j` of
ion species `i`, :math:`I_i^{(j)}` is the corresponding ionization rate,
:math:`R_i^{(j)}` the corresponding radiative recombination rate, and
:math:`\mathcal{I}_i^{(j)}` the corresponding fast electron impact ionization
coefficient. The difference between the "equilibrium" and "dynamic" options
above is merely in the transient term on the LHS. In equilibrium mode, the
ions are assumed to be in equilibrium, and so the transient term vanishes. In
the dynamic mode, on the other hand, transient effects are included.

Memory layout
-------------
As mentioned above, all ion species and charge states are stored in a single
unknown quantity called ``n_i``. Given :math:`N` ion species with atomic charges
:math:`Z_i\in\{ Z_1, Z_2, \ldots, Z_N \}`, defined on a radial grid with
:math:`n_r` grid points, ``n_i`` is stored as a single vector with
:math:`N_{n_i} = n_r\left[ \sum_i (Z_i+1) \right]` elements. The vector is laid
out so that :math:`n_r` adjacent elements in memory correspond to the same
radial density (i.e. to the same ion charge state).

A more graphical illustration is presented below.

+--------------+-------------------+----------------------------+-----------------------------+-----+
| **Address**  | 0-:math:`n_r-1`   | :math:`n_r`-:math:`2n_r-1` | :math:`2n_r`-:math:`3n_r-1` | ... |
+--------------+-------------------+----------------------------+-----------------------------+-----+
| **Contents** | :math:`n_1^{(0)}` | :math:`n_1^{(1)}`          | :math:`n_1^{(2)}`           | ... |
+--------------+-------------------+----------------------------+-----------------------------+-----+

The IonHandler class
--------------------
Since the unknown quantity ``n_i`` is stored in the same way as every any other
unknown quantity in DREAM, we cannot directly store information about the layout
of the ion densities with that unknown quantity. Instead, an ``IonHandler``
object is introduced for that purpose. The ``IonHandler`` keeps tracks of ion
names, atomic charges, and their order in the ``n_i`` solution vector. Via calls
to ``GetIndex(iZ, Z0)``, we can then translate a set consisting of an ion
species index (``iZ``) and the charge state number (``Z0``) to an index into the
``n_i`` vector containing the actual ion densities. The example below shows how
to access the radial density :math:`n_3^{(5)}`, corresponding to the fifth
charge state of the ion with index 3:

.. code-block:: c++

   const real_t *n_i = eqsys->GetUnknownData(id_n_i);
   len_t n35offs     = ionHandler->GetIndex(3, 5);

   const real_t *n35 = n_i + n35offs;

Ion equations
-------------
Due to the monolithic way in which ion densities are stored, regular equation
terms cannot be directly applied to ``n_i``. Instead, a class called
``IonEquationTerm`` has been derived to accomodate this need. The
``IonEquationTerm`` has all the same methods as the regular ``EquationTerm``,
but replaces the ``SetJacobianBlock()``, ``SetMatrixElements()`` and
``SetVectorElements()`` methods with the more specific
``SetCSJacobianBlock()``, ``SetCSMatrixElements()`` and
``SetCSVectorElements()``. The ``CS`` in the names of these methods stands for
*charge state*, and hence these methods set the matrix/vector elements for a
specified ion species charge state. The regular ``SetXXX()`` methods are as
usual called by the ``EquationSystem`` class when it is time to rebuild the
equation system matrix/vector, and these methods in turn iterate over all the
ion species charge states to which the equation term applies and calls the
``SetCSXXX()`` methods.

In addition to the regular input for the ``SetXXX()`` methods, the
``SetCSXXX()`` take three additional arguments as input:

- ``len_t iIon`` --- the ion species identifier
- ``len_t Z0`` --- the ion charge state number
- ``len_t rOffset`` --- the offset in the ``n_i`` vector that should be added in
  order to get the density corresponding to this ion.

When implementing a new ion equation term, it is the ``SetCSXXX()`` methods
which should be overridden, instead of the usual ``SetXXX()`` methods for the
regular ``EquationTerm``. Apart from this, the same rules apply to any
``IonEquationTerm`` as for the regular ``EquationTerm``, i.e. you will need to
also override the ``Rebuild()`` and ``GetNumberOfNonZerosPerRow()`` methods.

