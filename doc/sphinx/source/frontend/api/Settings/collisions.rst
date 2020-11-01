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
   everywhere for test-particles and supports the relativistic Maxwell-Jüttner 
   distribution as equilibrium, and is therefore the least approximate model
   currently implemented in DREAM. This mode is required for running
   fully-kinetic simulations where all electrons are resolved on the kinetic grid.

   The **super-thermal** limit of the full expressions has also been
   implemented and is intended for "standard" DREAM calculations, i.e. where
   a combined fluid-kinetic calculation is conducted. This model for the 
   distribution is non-conservative, where electrons are lost in the origin,
   are removed from the distribution and enter the fluid cold population.
   Notably, this setting removes the possibility of resolving Dreicer runaway
   generation, but captures all other mechanisms. 

   The **ultra-relativistic** is finally intended for benchmarking and for
   modelling the highly relativistic runaway electrons, where the :math:`v=c`
   limit is taken in the collision frequencies.

.. py:attribute:: collfreq_type

   Which type of ion nuclear screening to assume for ions.

   +---------------------------------------+-----------------------------------------------------------+
   | ``COLLFREQ_TYPE_COMPLETELY_SCREENED`` | Assume that ion charges are completely screened.          |
   +---------------------------------------+-----------------------------------------------------------+
   | ``COLLFREQ_TYPE_NON_SCREENED``        | Assume collisions with the full nuclear charge.           |
   +---------------------------------------+-----------------------------------------------------------+
   | ``COLLFREQ_TYPE_PARTIALLY_SCREENED``  | Take into account the partial screening of the nucleus    |
   +---------------------------------------+-----------------------------------------------------------+
   
   The **partially-screened** model is valid at all energies and captures
   the momentum-dependent screening effect of bound electrons, and contains the 
   most accurate collision model available. Formally,
   the **completely-screened** type is the low-energy limit of this general
   model, and **non-screened** the asymptotic high-energy limit.

.. py:attribute:: lnlambda

   Which model to use for the Coulomb logarithm :math:`\ln\Lambda`.

   +---------------------------------------+-----------------------------------------------------------------------+
   | ``LNLAMBDA_THERMAL``                  | Take the Coulomb logarithm to depend only on temperature and density. |
   +---------------------------------------+-----------------------------------------------------------------------+
   | ``LNLAMBDA_CONSTANT``                 | Take the Coulomb logarithm to depend only on temperature and density. |
   +---------------------------------------+-----------------------------------------------------------------------+
   | ``LNLAMBDA_ENERGY_DEPENDENT``         | Account for the Coulomb logarithm's dependence on electron energy.    |
   +---------------------------------------+-----------------------------------------------------------------------+

   The energy-dependent Coulomb logarithm accounts for the increase of (and difference between) 
   the electron-electron and electron-ion Coulomb logarithms as the electron energy increases.
   ``LNLAMBDA_THERMAL`` corresponds roughly to the energy-dependent model evaluated at the 
   thermal momentum, and ``LNLAMBDA_CONSTANT`` at the speed of light.

.. py:attribute:: pstar_mode

   Which model to use for the effective critical momentum :math:`p_\star`, which set 
   the analytical runaway rates of avalanche, Compton and tritium.

   +--------------------------------+----------------------------------------------------------------------+
   | ``PSTAR_MODE_COLLISIONLESS``   |  Assumes banana-regime runaway generation with trapping corrections  |
   +--------------------------------+----------------------------------------------------------------------+
   | ``PSTAR_MODE_COLLISIONAL``     |  Assumes local runaway generation                                    |
   +--------------------------------+----------------------------------------------------------------------+

   The local model is inspired by C J McDevitt, X-Z Tang, EPL 127, 45001 (2019), where 
   they showed that at sufficiently strong electric fields (and therefore low critical momenta),
   the electrons enter the collisional regimes where they are scattered before completing
   a bounce orbit. At electron densities of 2e20, this occurs at approximately :math:`p_c\sim 0.1`. 

   .. note::
      Fluid Dreicer generation does not support `COLLISIONLESS` and always default to `COLLISIONAL`.

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

