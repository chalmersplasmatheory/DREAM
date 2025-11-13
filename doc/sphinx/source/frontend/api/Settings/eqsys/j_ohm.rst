.. _ds-eqsys-j_ohm:

Ohmic current density
=====================
The ohmic current density is tracked with the unknown ``j_ohm``, and is always
calculated self-consistently from other unknowns.

Unless the ohmic current is resolved on the kinetic grid, it is given by

.. math::
    \frac{j_\mathrm{ohm}}{B} 
    = \sigma \frac{ \langle \boldsymbol{E}\cdot\boldsymbol{B} \rangle }{  \langle B^2 \rangle  }

where :math:`\sigma` denotes the conductivity, and depends on the plasma parameters 
(most sensitively to the cold electron temperature ``T_cold``). When the ohmic current
is resolved by the distribution function, it is defined as

.. math::
    \frac{j_\mathrm{ohm} + j_\mathrm{hot}}{B} = \frac{2\pi e}{B_\mathrm{min}} 
    \int \mathrm{d}p \mathrm{d}\xi_0 \,\frac{\mathcal{V}'}{V'} \Theta v f_\mathrm{hot},

where :math:`\Theta` is zero for trapped particles and equals unity for passing particles,
and the contribution from the hot-electron current is denoted ``j_hot`` (typically defined as 
a similar integral moment, but integrated for momenta above some cutoff). 

.. contents:: Page overview
   :local:
   :depth: 3



.. _ds-eqsys-j_ohm-conductivity:

Conductivity formulas
---------------------

The conductivity can be calculated using different models for the conductivity. 
DREAM implements the model described in 
`O Sauter, C Angioni and YR Lin-Liu (1999) <https://doi.org/10.1063/1.873240>`_,
which is a numerical fit of the neoclassical conductivity to Fokker-Planck simulations
across all collisionality regimes (and not just the banana regime, which is available
with kinetic simulations in DREAM). In the Sauter paper, the classical conductivity 
(in the absence of neoclassical effects, i.e. in cylindrical geometry) appears as a free parameter, 
which we have chosen as the relativistic conductivity given in 
`BJ Braams and CFF Karney (1989) <https://doi.org/10.1063/1.858966>`_,
where we interpolate in the values they provide in Table 1.

The conductivity model can be controlled with the following settings:

+----------------------------------------------+-------------------------------------------------------------------------------------------------+
| Name                                         | Description                                                                                     |
+==============================================+=================================================================================================+
| ``CONDUCTIVITY_MODE_BRAAMS``                 | Uses the cylindrical conductivity, ignoring all neoclassical effects                            |
+----------------------------------------------+-------------------------------------------------------------------------------------------------+
| ``CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS``   | Uses the Sauter formula in the collisionless limit, representing the banana-regime conductivity |
+----------------------------------------------+-------------------------------------------------------------------------------------------------+
| ``CONDUCTIVITY_MODE_SAUTER_COLLISIONAL``     | Uses the full expression given in Eqs (13a-b) in the Sauter paper                               |
+----------------------------------------------+-------------------------------------------------------------------------------------------------+

.. note::
    The default conductivity mode in DREAM is ``CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS`` 
    as this corresponds to the same approximations as used in the kinetic equation.  


Corrected conductivity
----------------------

When resolving ohmic current kinetically in DREAM, an order-unity error is incurred 
from the fact that we do not include the electron-electron field-particle collision operator.
At an effective charge of :math:`Z_\mathrm{eff} = 1`, the error is approximately a factor of 2,
but approaches 1 as the effective charge goes to infinity.

In order to correct for this error, DREAM allows the option to add an ohmic current correction
on top of the kinetic contribution to ``j_ohm``. An extensive numerical scan of kinetic 
simulations has revealed that the ratio of numerical conductivity (obtained in DREAM with a 
test-particle collision operator) to  **real** (as described by the Braams-Karney model) is 
well described by

.. math::
    \frac{\sigma_\mathrm{numerical}}{\sigma} \approx 1 + \frac{1}{b+Z_\mathrm{eff}},

where we numerically determined :math:`a=-1.406` and :math:`b=1.888`, yielding a 
maximum relative error less than :math:`1\%` across a large range of plasma temperatures
and charge.

The conductivity corretion used in DREAM is controlled with the following settings:

+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Name                                | Description                                                                                                                                                                             |
+=====================================+=========================================================================================================================================================================================+
| ``CORRECTED_CONDUCTIVITY_DISABLED`` | No correction is added.                                                                                                                                                                 |
+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``CORRECTED_CONDUCTIVITY_ENABLED``  | A correction is added of the form :math:`\Delta j_\mathrm{hot} = (\sigma - \sigma_\mathrm{numerical})\frac{ \langle \boldsymbol{E}\cdot\boldsymbol{B} \rangle }{  \langle B^2 \rangle}` |
+-------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


.. note::
    The conductivity correction is only active when the hot-tail grid is enabled,
    when using ``COLLISION_FREQENCY_MODE_FULL`` and when resolving the pitch-angle distribution. 



Example
*******

The conductivity mode and correction can be set according to the following example:

.. code-block:: python
    
    import DREAM.Settings.Equations.OhmicCurrent as JOhm

    ds = DREAMSettings()

    ds.eqsys.j_ohm.setConductivityMode(JOhm.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)
    ds.eqsys.j_ohm.setCorrectedConductivity(JOhm.CORRECTED_CONDUCTIVITY_ENABLED)


Set initial current density
---------------------------
.. _ds-eqsys-j_ohm-init:
To solve Faraday's law, an initial condition on :math:`E_\parallel` is required.
This is usually taken to be a value such that a desired current density profile
is obtained. Since this is such a common way of initializing the electric field,
functionality exists in DREAM to initialize :math:`E_\parallel` based on a given
current density profile.

To initialize the electric field :math:`E_\parallel` according to a given
current density profile, you can use the ``j_ohm.setInitialProfile()`` method:

.. code-block:: python

    import numpy as np

    ds = DREAMSettings()
    ...

    # Plasma minor radius
    a = 1.5
    r = np.linspace(0, a)
    j0 = 1e6 * (1 - (1-0.001**(1/0.41))*(r/a)**2)**0.41

    # Set initial value for j_tot
    ds.eqsys.j_ohm.setInitialProfile(j0, radius=r)

.. warning::

   Although the unknown quantity which the setting is made on is ``j_ohm``, it
   is actually the total current density profile ``j_tot`` which is set with
   the above code. This means that any initial runaway current is also counted
   as part of the value set with ``j_ohm.setInitialProfile(j0 ,r)``. The actual
   ohmic current density will thus be ``j0 - j_re``.

Re-scaling profile to Ip
************************
Another common situation is that in which you know the shape of the desired
current profile, and what total current :math:`I_{\rm p}` it should integrate
to, but not the appropriate normalization coefficient to use to achieve this.
In this situation, the ``j_ohm.setInitialProfile()`` function provides the
``Ip0`` argument, which allows you to specify the plasma current to which the
profile should integrate. DREAM will then calculate the appropriate
normalization to use in the kernel, and normalize the profile accordingly:

.. code-block:: python
    
    ...
    Ip0 = 10e6  # 10 MA current
    ds.eqsys.j_ohm.setInitialProfile(j0, r, Ip0=Ip0)

Initializing from other current densities
*****************************************
Current densities in DREAM are formally the parallel current density in the
point where :math:`B=B_{\rm min}` along a flux surface. Other codes, or
experiments, tend to work with different current densities, and as such one must
be cautious when using current density profiles obtained from other
codes/experiments in DREAM. Since other definitions of current densities are
usually related to the DREAM definition via a product of geometrical factors,
DREAM can automatically re-scale the given current density profile, assuming it
uses one of the supported definitions.

The definitions of current density supported (for the initial current density
profile) in DREAM are as follows:

+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------+
| **Option name**                | **Mathematical expression**                                                                                                                                      | **Comment**                                        |
+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------+
| ``PROFILE_TYPE_J_PARALLEL``    | :math:`j_\parallel = j_\parallel`                                                                                                                                | Parallel current density. Same as used internally. |
+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------+
| ``PROFILE_TYPE_J_DOT_GRADPHI`` | :math:`\left\langle\boldsymbol{j}\cdot\nabla\phi\right\rangle = \frac{j_{\parallel,\rm min}}{B_{\rm min}}\left[ G\left\langle\frac{1}{R^2}\right\rangle \right]` | Flux-surface averaged toroidal current density.    |
+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------+
| ``PROFILE_TYPE_JTOR_OVER_R``   | :math:`\bar{J}_\phi = \frac{j_{\parallel,\rm min}}{B_{\rm min}}\frac{\left\langle B^2\right\rangle}{G\left\langle 1/R^2\right\rangle}`                           | Normalized toroidal current density.               |
+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------+
| ``PROFILE_TYPE_CORSICA``       | (see above)                                                                                                                                                      | Synonym for ``PROFILE_TYPE_JTOR_OVER_R``.          |
+--------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------+

To use these re-scalings, just add the argument ``profile_type`` with one of the
options above:

.. code-block:: python

    import DREAM.Settings.Equations.OhmicCurrent as JOhm
    ...
    ds.eqsys.j_ohm.setInitialProfile(
        j0, r, Ip0,
        profile_type=JOhm.PROFILE_TYPE_CORSICA
    )


