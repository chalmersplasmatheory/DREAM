.. _ds-collisions:

collisions
==========
The ``collisions`` property of the ``DREAMSettings`` object contains settings
for how electron-electron and electron-ion collisions are to be handled. The
object is very simple in structure and only contains a number of flags. The
available options are listed below.

Object documentation
--------------------
.. py:class:: CollisionHandler

.. py:attribute:: bremsstrahlung_mode

   Whether or not to take bremsstrahlung losses into account. The available
   options are

   +----------------------------------------+---------------------------------+
   | ``BREMSSTRAHLUNG_MODE_NEGLECT``        | Neglect bremsstrahlung losses   |
   +----------------------------------------+---------------------------------+
   | ``BREMSSTRAHLUNG_MODE_STOPPING_POWER`` | Use the stopping-power formula. |
   +----------------------------------------+---------------------------------+

.. py:attribute:: collfreq_mode

   Which model to use for the collision frequencies. The available options are

   +--------------------------------------+----------------------------------------------------+
   | ``COLLFREQ_MODE_SUPERTHERMAL``       | Super-thermal limit of collision frequencies.      |
   +--------------------------------------+----------------------------------------------------+
   | ``COLLFREQ_MODE_FULL``               | Fully relativistic collision frequencies.          |
   +--------------------------------------+----------------------------------------------------+
   | ``COLLFREQ_MODE_ULTRA_RELATIVISTIC`` | Ultra-relativistic limit of collision frequencies. |
   +--------------------------------------+----------------------------------------------------+

   The **full** collision model uses collision frequencies which are valid 
   everywhere for test-particles, and is therefore the least approximate model
   currently implemented in DREAM. This mode is required for running
   fully-kinetic simulations where all electrons are kinetic (and the cold
   electron density ``n_cold`` is identically zero).

   The **super-thermal** limit of the full expressions has also been
   implemented and is intended for "standard" DREAM calculations, i.e. where
   a combined fluid-kinetic calculation is conducted. The super-thermal
   collision frequencies are then applied to the hot electrons, which are
   evolved kinetically at energies of a several tens of keV.

   The **ultra-relativistic** is finally intended for benchmarking and for
   modelling the highly relativistic runaway electrons.

.. py:attribute:: collfreq_type

   Which type of ion nuclear screening to assume for ions.

   +---------------------------------------+-----------------------------------------------------------+
   | ``COLLFREQ_TYPE_COMPLETELY_SCREENED`` | Assume that ion charges are completely screened.          |
   +---------------------------------------+-----------------------------------------------------------+
   | ``COLLFREQ_TYPE_NON_SCREENED``        | Assume that full ion charge can be probed.                |
   +---------------------------------------+-----------------------------------------------------------+
   | ``COLLFREQ_TYPE_PARTIALLY_SCREENED``  | Take into account the partial screening of the ion charge |
   +---------------------------------------+-----------------------------------------------------------+

   The **completely-screened** limit is valid in weakly-ionized (e.g.
   low-temperature) plasmas, while the **non-screened** limit applies to fully
   ionized (e.g. high-temperature) plasmas. The **partially-screened** model is
   valid in both these limits, as well as between them, and is therefore the
   most complete screening model.

.. py:attribute:: lnlambda

   Which model to use for the Coulomb logarithm :math:`\ln\Lambda`.

   +---------------------------------------+-----------------------------------------------------------------------+
   | ``LNLAMBDA_CONSTANT``                 | Take the Coulomb logarithm to depend only on temperature and density. |
   +---------------------------------------+-----------------------------------------------------------------------+
   | ``LNLAMBDA_ENERGY_DEPENDENT``         | Account for the Coulomb logarithm's dependence on electron energy.    |
   +---------------------------------------+-----------------------------------------------------------------------+

   The energy-dependent Coulomb logarithm becomes important in partially and
   weakly ionized plasmas and when the electrons become highly energetic.

.. py:attribute:: pstar_mode

   Which model to use for the effective critical momentum :math:`p_\star`.

   +---------------------------------------+-----------------------------------------------------------------------+
   | ``PSTAR_MODE_COLLISIONAL``            |                                                                       |
   +---------------------------------------+-----------------------------------------------------------------------+
   | ``PSTAR_MODE_COLLISIONLESS``          |                                                                       |
   +---------------------------------------+-----------------------------------------------------------------------+

.. todo::

   Describe the options for ``pstar_mode``.

Examples 
--------
The below example shows how to set collision options:

.. code-block:: python

   from DREAM.DREAMSettings import DREAMSettings
   import DREAM.Settings.CollisionHandler as Collisions

   ds = DREAMSettings()
   ...
   ds.collisions.bremsstrahlung_mode = Collisions.BREMSSTRAHLUNG_MODE_STOPPING_POWER
   ds.collisions.collfreq_mode       = Collisions.COLLFREQ_MODE_FULL
   ds.collisions.collfreq_type       = Collisions.COLLFREQ_TYPE_NON_SCREENED
   ds.collisions.lnlambda            = Collisions.LNLAMBDA_CONSTANT
   ds.collision.pstar_mode           = Collisions.PSTAR_MODE_COLLISIONAL

