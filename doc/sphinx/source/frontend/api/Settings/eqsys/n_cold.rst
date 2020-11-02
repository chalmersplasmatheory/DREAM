.. _ds-eqsys-n_cold:

Cold electron density
=====================
The density of (free) cold electrons is tracked with the unknown ``n_cold``.
This parameter can either be prescribed or calculated self-consistently.

.. note::

   Unless explicitly specified, the cold electron density is evolved
   self-consistently.

.. contents:: Page overview
   :local:
   :depth: 3

Prescribed evolution
--------------------
The cold electron density can be prescribed in time and space using the method
``setPrescribedData()``.

.. todo::

   Should we disable prescribing ``n_cold`` when the hot electron grid is
   enabled?

Example
*******
The following example illustrates how to prescribe the cold electron density
in a DREAM simulation:

.. code-block:: python

   ds = DREAMSettings()
   ...
   r = np.linspace(r0, r1, nr)
   t = np.linspace(t0, t1, nt)

   # Generate density profiles
   n = np.array([...])      # shape (nt, nr)

   # Prescribe cold electron density
   ds.eqsys.n_cold.setPrescribedData(density=n, radius=r, times=t)

A constant and uniform density can also be prescribed with the more compact
syntax

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.eqsys.n_cold.setPrescribedData(density=2e19)


Self-consistent evolution
-------------------------
When evolved self-consistently, the cold electron density is obtained from the
electron conservation equation

.. math::

   n_{\rm cold} = n_{\rm free} - n_{\rm hot} - n_{\rm RE},

where :math:`n_{\rm free}` is the density of free electrons, :math:`n_{\rm hot}`
is the density of electrons on the hot electron grid (if enabled) and
:math:`n_{\rm RE}` is the density of runaway electrons. To evolve the cold
electrons self-consistently, use the option
``ColdElectrons.TYPE_SELFCONSISTENT``.

Example
*******
The following example illustrates how to evolve the cold electron density
self-consistently in a DREAM simulation:

.. code-block:: python

   import DREAM.Settings.Equations.ColdElectrons as ColdElectrons

   ds = DREAMSettings()
   ...
   ds.eqsys.n_cold.setType(ColdElectrons.TYPE_SELFCONSISTENT)


Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.ColdElectrons.ColdElectrons
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

