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
   possible to set up a simulation in which ions and electrons are unbalanced,
   **thus breaking quasi-neutrality**.
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
+====================================+=================================================================================================================================+
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


Adding ions
-----------
Ion species are added to a DREAM simulation one by one using the ``addIon()``
method. When adding a new ion species, one must specify the name of the species
(which must be unique), the atomic charge :math:`Z`, which equation to use for
evolving the ion densities, as well as any initial or prescribed ion densities.
Since DREAM uses a user-specified name parameter to keep track of ion species,
it is possible to add several ion populations of the same species and evolve
them independently of each other.

When initializing ion densities it is also possible to specify which charge
state :math:`Z_0` should have a non-zero density. This is done using the
``adIon()`` input parameter ``Z0``. If ``Z0`` is not specified, and the
``iontype`` parameter does not specify the charge state, the user must provide
the input density ``n`` as an array with each element of the first dimension
representing the density in the corresponding charge state.

When adding a tritium ion species the ``tritium`` parameter should also be set
to ``True``. This allows DREAM to correctly identify the species as tritium and
use its density when calculating the
:ref:`tritium decay runaway rate<ds-eqsys-ions-tritium>`.

.. warning::

   If you would like to add ions to a :ref:`restart<restart>` simulation then
   you will not be able to initialize the ion densities from the previous
   output. Instead you will have to prescribe the initial ion densities again.

   *The reason why ion densities cannot be initialized from a previous DREAM
   output if new ions are added is because DREAM is not able to determine if
   ions are added, removed or renamed. This makes it impossible for DREAM to
   keep track of which ion densities in the previous output map to which ions in
   the new input.*

Example
^^^^^^^
Prescribed ions
...............
The following example illustrates how to prescribe the density to be constant
and uniform for an ion species:

.. code-block:: python

   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ...
   ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_PRESCRIBED, Z0=1, n=2e19)

For more advanced use cases, one can also prescribe the ion densities for each
charge state, at a set of radii and at a set of time points.

.. code-block:: python

   import numpy as np
   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ...
   # Ion name and atomic charge
   name = 'Ar'
   Z    = 18

   # Define radial and time grids
   r = np.linspace(r0, r1, nr)
   t = np.linspace(t0, t1, nt)

   # Construct ion densities
   n = np.array([...])      # Shape (Z+1, nt, nr)

   # n[0,:] = Neutral charge state
   # ...
   # n[Z,:] = Fully ionized charge state

   ds.eqsys.n_i.addIon(name=name, Z=Z, iontype=Ions.IONS_PRESCRIBED, n=n, r=r, t=t)

.. note::

   All ion densities must be prescribed on the same radial/time grid. Scalar
   ions densities can still be prescribed, must be so *after* prescribing at
   least one ion species with the non-scalar radial/time grid.

Dynamic ions
............
Dynamic ions can be added similarly to prescribed ions, but they are instead
evolved using the ion rate equations. The following example illustrates how to
add an ion species that is evolved using the ion rate equations and is
initialized with a uniform radial density profile with all ions in the fully
ionized charge state:

.. code-block:: python

   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ...
   ds.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, Z0=1, n=2e19)

Alternatively, the ion charge state densities can be directly prescribed. The
following example illustrates how to prescribe the initial ion charge state
density profiles and evolve them using the ion rate equations:

.. code-block:: python

   import numpy as np
   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ...
   # Ion name and atomic charge
   name = 'Ar'
   Z    = 18

   # Define radial grid
   r = np.linspace(r0, r1, nr)

   # Construct ion densities
   n = np.array([...])      # Shape (Z+1, nr)

   # n[0,:] = Neutral charge state
   # ...
   # n[Z,:] = Fully ionized charge state

   ds.eqsys.n_i.addIon(name=name, Z=Z, iontype=Ions.IONS_DYNAMIC, n=n, r=r)

.. note::

   All ion densities must be prescribed on the same radial/time grid. Scalar
   ions densities can still be prescribed, must be so *after* prescribing at
   least one ion species with the non-scalar radial/time grid.



Ion temperature
---------------
When the temperature in the plasma is solved for self-consistently, the default 
mode is to model only the electron temperature evolution due to radiation, heating
and transport. It is possible to also include the evolution of the temperature of 
each ion species, where different charge states of the same species are assumed to 
have the same temperature. Ions and electrons exchange energy via elastic collisions,
taking the form of a rate equation

.. math::

   \frac{\partial W_i}{\partial t} = \sum_j Q_{ij} + Q_{ie}, 

where the sum over :math:`j` is taken over all ion species, and :math:`Q_{ij}` 
is the collisional energy transfer integrated over Maxwellians of different 
density :math:`N_i = \sum_j n_i^{(j)}` and heat :math:`W_i = 3 e N_i T_i / 2`,

.. math::

   Q_{ij} = \sqrt{\frac{3}{\pi}} \frac{Z_i^2 Z_j^2 e^4 \ln\Lambda_{ij} N_i N_j\sqrt{m_i N_i m_j N_j} }{4\pi \varepsilon_0^2}
   \frac{N_i W_j - N_j W_i}{(m_j N_j W_i + m_i N_i W_j)^{3/2}}.

In addition, the cold-electron temperature will pick up a similar contribution,

.. math::

   \left( \frac{\partial W_\mathrm{cold}}{\partial t} \right)_Q = \sum_i Q_{ei}

where the anti-symmetry of :math:`Q_{ij} = -Q_{ji}` ensures that the sum of the heat equations
for all ions and electrons exhibits energy conservation by these elastic collisions.

The initialization and behavior of ion temperatures are controlled via the ``T`` argument
in ``Ions::addIon(..., T)``. It behaves as follows:

- If at least one ``T`` is explicitly set, DREAM will add the additional quantities ``N_i`` and ``W_i`` to the equation system and evolve them as ``nontrivial`` unknowns.
- Species which are not explicitly set will be initialized to ``T=0``.
- If the ``type`` for ``T_cold`` is set to ``TYPE_SELFCONSISTENT``, the ion and electron heat ``W_i`` and ``W_cold`` will be evolved according to the equations above. If ``TYPE_PRESCRIBED`` it will be given by its initial value ``W_i=constant``.


Example
^^^^^^^
In the below example, we consider a scenario where a hydrogenic population from before
a disruption has the same initial temperature as the electrons, whereas an injected 
hydrogenic species and Neon impurity is introduced at effectively zero temperature.
The ions and electrons will be allowed to exchange energy via collisions.
All ion and electron initial temperatures are taken to be uniform in radius.

.. code-block:: python

   import numpy as np
   import DREAM.Settings.Equations.IonSpecies as Ions
   import DREAM.Settings.Equations.ColdElectronTemperature as T_cold

   ds = DREAMSettings()

   ...

   T_initial = 4e3 #eV

   ds.eqsys.T_cold.setPrescribedData(T_initial)
   ds.eqsys.n_i.addIon(name='D',     Z=1,  T=T_initial, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED, n=1e20)
   ds.eqsys.n_i.addIon(name='Ne',    Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e19)
   ds.eqsys.n_i.addIon(name='D_inj', Z=1,  iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e21)

   ds.eqsys.T_cold.setType(ttype=T_cold.TYPE_SELFCONSISTENT)


Atomic data
-----------
Various types of atomic data are downloaded from ADAS and NIST during the
configuration phase of DREAM (prior to compilation). The data is downloaded
using the scripts ``tools/get_adas.py`` and ``tools/get_nist.py``, which should
be automatically invoked by CMake.

ADAS
^^^^
  *The Atomic Data and Analysis Structure (ADAS) is an interconnected set of
  computer codes and data collections for modelling the radiating properties of
  ions and atoms in plasmas. It can address plasmas ranging from the 
  interstellar medium through the solar atmosphere and laboratory thermonuclear
  fusion devices to technological plasmas. ADAS assists in the analysis and
  interpretation of spectral emission and supports detailed plasma models.*

  -- ADAS website, https://www.adas.ac.uk/about.php

DREAM utilizes data for four types of coefficients from the OpenADAS database
(https://open.adas.ac.uk), name ACD (effective recombination coefficients),
SCD (effective ionization coefficients), PLT (line power driven by excitation of
dominant ions) and PRB (continuum and line power driven by recombination and
bremsstrahlung of dominant ions). Data is downloaded for the elements specified
in the file ``tools/elements.json``. The configuration file is given in JSON
format with one entry per element in a key-value format. The key should be the
official name of the element and the value is the year in which the dataset to
use was published.

The ``get_adas.py`` script can also be run manually, separately from the DREAM
CMake script to generate a C/C++ file of ADAS data. The script is then invoked
from the command line and accepts the following arguments:

+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| Name                 | Description                                                                                                                 |
+======================+=============================================================================================================================+
| ``--cachedir DIR``   | Use ``DIR`` to store raw data downloaded from ADAS. This can be re-used later to avoid making an HTTP request to Open-ADAS. |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--elements FILE``  | Load elements to fetch data for from the JSON file ``FILE``.                                                                |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--no-cache``       | Do **not** store raw ADAS data locally.                                                                                     |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--no-compile``     | Do **not** generate C++ source files with the rate coefficients.                                                            |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``-o``, ``--output`` | Name of output C++ source file to generate.                                                                                 |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--type-int``       | C++ type to use for integers.                                                                                               |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--type-real``      | C++ type to use for real numbers.                                                                                           |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+

NIST
^^^^
DREAM uses ionization and binding energy values tabulated in the database of the
american National Institute of Standards (NIST) database. The values are
downloaded by making a request via the form
https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html and are then embedded
into a C++ source file. The elements to download are specifed by name in an
array in the ``main()`` function of the script ``tools/get_nist.py``.

While the ``tools/get_nist.py`` script should be automatically invoked when
running CMake for DREAM, it can also be run manually. The script accepts the
following command-line arguments:

+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| Name                 | Description                                                                                                                 |
+======================+=============================================================================================================================+
| ``--cachedir DIR``   | Use ``DIR`` to store raw data downloaded from ADAS. This can be re-used later to avoid making an HTTP request to Open-ADAS. |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--ionization``     | Downloads ionization energy data instead of the default binding energy data.                                                |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--no-cache``       | Do **not** store raw ADAS data locally.                                                                                     |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--no-compile``     | Do **not** generate C++ source files with the rate coefficients.                                                            |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``-o``, ``--output`` | Name of output C++ source file to generate.                                                                                 |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--type-int``       | C++ type to use for integers.                                                                                               |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+
| ``--type-real``      | C++ type to use for real numbers.                                                                                           |
+----------------------+-----------------------------------------------------------------------------------------------------------------------------+


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


.. _ds-eqsys-ions-tritium:

Tritium
-------
Runaway electrons can be generated when tritium decays into helium-3 via the
beta decay process

.. math::

   \mathrm{T} \to \ ^3_2\mathrm{He} + \mathrm{e}^- + \bar{\nu}_{\rm e}.

The corresponding runaway rate is given by

.. math::

   \left( \frac{\mathrm{d} n_{\rm RE}}{\mathrm{d} t} \right)_{\rm T} \approx
   \ln 2 \frac{n_{\rm T}}{\tau_{\rm T}} F_\beta\left( \gamma_{\rm c} \right),

where :math:`n_{\rm T}` is the tritium density, :math:`\tau_{\rm T} = 4800\pm 8`
days is the tritium half-life, and :math:`F_\beta(\gamma_{\rm c})` denotes the
fraction of beta electrons generated with an energy above the critical energy
:math:`\gamma_{\rm c}` for runaway to occur.

Tritium runaway generation is enabled with the :ref:`n_re<ds-eqsys-n_re>`
object, but it is necessary to also provide a tritium ion species in the
ions interface. The user must make sure to specify ``tritium=True`` when adding
the tritium ion to the simulation using the ``addIon()`` method, as described
above under `Adding ions`_.

.. note::

   It is possible to include multiple tritium populations in the simulation,
   simply by providing ``tritium=True`` when adding each of them.

Example
^^^^^^^
The following example illustrates how to enable the tritium decay runaway
mechanism in a DREAM simulation:

.. code-block:: python

   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ...
   # Include source term in equation for n_re
   ds.eqsys.n_re.setTritium(True)

   # Add tritium ion species to list of ions
   ds.eqsys.n_i.addIon('T', Z=1, iontype=Ions.IONS_DYNAMIC, n=2e19, tritium=True)


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

