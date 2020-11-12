.. _ds-eqsys-n_re:

RunawayElectrons
================
The ``RunawayElectons`` object, located under ``.eqsys.n_re`` in the
:ref:`DREAMSettings` object is used to specify settings for the fluid quantity
``n_re``, representing the density of runaway electrons in simulations. The
quantity is always present in DREAM simulations and is evolved using the
equation

.. math::

   \frac{\mathrm{d} n_{\rm RE}}{\mathrm{d} t} =
   \gamma_{\rm Dreicer} + \Gamma_{\rm Ava}n_{\rm RE} + \Phi^{(p)}_{\rm hot} +
   \gamma_{\rm tritium} + \gamma_{\rm Compton} +
   \frac{1}{V'}\frac{\partial}{\partial r}\left[
       V'\left( An_{\rm RE} +  D\frac{\partial n_{\rm RE}}{\partial r} \right)
   \right]

where

- :math:`\gamma_{\rm Dreicer}` is the runaway rate due to Dreicer mechanism.
- :math:`\Gamma_{\rm Ava}` is the avalanche growth rate.
- :math:`\Phi^{(p)}_{\rm hot}` is the flux of particles from the (kinetic) hot-tail grid.
- :math:`\gamma_{\rm tritium}` is the runaway rate due to tritium decay.
- :math:`\gamma_{\rm Compton}` is the runaway rate due to Compton scattering.

and the last term represents advective-diffusive radial transport. Each term can
be individually enabled/disabled using the methods described further down on
this page.

.. note::

   The Dreicer generation rate :math:`\gamma_{\rm Dreicer}` should typically
   be disabled in simulations including the hot-tail grid. This is because the
   Dreicer generation will then be modelled using a kinetic equation and the
   Dreicer generation rate become a part of the hot electron flux
   :math:`\Phi^{(p)}_{\rm hot}`.


.. note::

   Note that the runaway density :math:`n_{\rm RE}` is **never** defined in
   terms of the runaway electron distribution function :math:`f_{\rm RE}` in
   DREAM. The quantities are instead evolved separately and set up in such a
   way that their local densities should agree in the end if the simulation was
   set up in a consistent way.

   *It is possible---and perfectly okay---for* :math:`n_{\rm RE}` *to differ
   from the density moment of* :math:`f_{\rm RE}` *if the latter loses a
   significant number of particles at* :math:`p=p_{\rm max}`.


(Effective) critical electric field
-----------------------------------

.. todo::

   Describe the effective critical electric field settings.


Initialization
--------------
By default, the runaway density is set to zero at the beginning of a DREAM
simulation. It is however possible to start the simulation with a non-zero
runaway density profile, and this is done using the ``setInitialProfile()``
method.

.. warning::

   If the runaway grid is enabled, one should also make sure to appropriately
   initialize the runaway electron distribution function when setting a non-zero
   initial runaway density. The density moment of the runaway electron
   distribution function :math:`f_{\rm RE}` should always agree with the runaway
   electron density :math:`n_{\rm RE}`.

Example
^^^^^^^
The following example illustrates how to initialize the runaway electron
density:

.. code-block:: python

   import numpy as np

   ds = DREAMSettings()
   ...
   # Define radial grid for the initial profile
   r = np.linspace(r0, r1, nr)

   # Generate the runaway electron density profile
   nre = np.array([...])        # shape (nr,)

   ds.eqsys.n_re.setInitialProfile(density=nre, r=r)

If a uniform initial density profile is desired one can specify the ``density``
parameter as a scalar value instead:

.. code-block:: python

   import numpy as np

   ds = DREAMSettings()
   ...
   ds.eqsys.n_re.setInitialProfile(density=2e16)


Runaway generation rates
------------------------
DREAM includes a number of common runaway electron source terms for use in both
pure fluid as well as combined fluid-kinetic simulations.

Avalanche
^^^^^^^^^

.. todo::

   Write about avalanche generation in DREAM.


Compton
^^^^^^^

.. todo::

   Write about Compton runaway generation in DREAM.


Dreicer
^^^^^^^

.. todo::

   Write about Dreicer runaway generation in DREAM.


Tritium
^^^^^^^
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

Tritium runaway generation is enabled with the ``setTritium()`` method of the
``RunawayElectrons`` object, but it is necessary to also provide a tritium ion
species in the :ref:`ions interface <ds-eqsys-ions-tritium>`. The user must make
sure to specify ``tritium=True`` when adding the tritium ion to the simulation
using the ``addIon()`` method.

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
   ds.eqsys.ions.addIon('T', Z=1, iontype=Ions.IONS_DYNAMIC, n=2e19, tritium=True)




Class documentation
-------------------
.. autoclass:: DREAM.Settings.Equations.RunawayElectrons.RunawayElectrons
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__


