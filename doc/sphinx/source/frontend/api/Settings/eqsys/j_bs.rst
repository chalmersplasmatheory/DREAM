.. _ds-eqsys-j_bs:

Bootstrap current density
=====================
The bootstrap current density is tracked with the unknown ``j_bs``, and is disabled 
by default as it typically vanishes during disruptions as the plasma profiles flattens.
It is determined mainly by gradients in the the density and temperature profiles,
following the formula

.. math::
    \frac{j_\mathrm{bs}}{B}
    = -\frac{G}{ \langle B^2\rangle \psi'_\mathrm{ref} }
    \bigg[ \mathcal{L}_{31} \frac{p}{n} \frac{ \partial \ln n }{ \partial r }
    +( \mathcal{L}_{31} + \mathcal{L}_{32} ) p_e \frac{ \partial \ln T_e }{ \partial r } 
    +\mathcal{L}_{31} (1+\alpha) \sum_i p_i \frac{ \partial \ln T_i }{ \partial r } \bigg],

where :math:`p` and :math:`n` are the total pressure and density, respectively, 
summed over the thermal electrons and all ion species. In case ion temperatures 
are not evolved in time, then ions are assumed to be in thermal equilibrium with 
the electrons, and :math:`T_i=T_e`.
The coefficients :math:`\mathcal{L}_{31}`, :math:`\mathcal{L}_{32}`, and :math:`\alpha` are
functions of the local fraction of trapped particles, the effective charge, and
the electron collisionality (or ion collisionality for :math:`\alpha`).


.. contents:: Page overview
   :local:
   :depth: 3




Bootstrap models
----------------
.. _ds-eqsys-bs-models:

Various models exist, and these mainly vary in how the coefficients :math:`\mathcal{L}_{31}`, 
:math:`\mathcal{L}_{32}`, and :math:`\alpha` are defined. Currently, only the model given in
`Redl et al. (2021) <https://doi.org/10.1063/5.0012664>`_ is implemented in DREAM. 

The following table summarises the available options:

+----------------------------------------------+-------------------------------------------------------------------------------------------------+
| Name                                         | Description                                                                                     |
+==============================================+=================================================================================================+
| ``BOOTSTRAP_MODE_DISABLED``                  | Neglects the bootstrap current, setting :math:`j_\mathrm{bs}=0`.                                |
+----------------------------------------------+-------------------------------------------------------------------------------------------------+
| ``BOOTSTRAP_MODE_REDL``                      | Uses the model given in the Redl paper.                                                         |
+----------------------------------------------+-------------------------------------------------------------------------------------------------+


The bootstrap current model can be set following the following example:

.. code-block:: python

    import DREAM.Settings.Equations.BootstrapCurrent as BootstrapCurrent

    ds = DREAMSettings()

    ds.eqsys.j_bs.setMode(BootstrapCurrent.BOOTSTRAP_MODE_REDL)




Initialisation mode
-------------------
.. _ds-eqsys-bs-init:

When providing the initial current density in a simulation, it can either be the
total current density :math:`j_\mathrm{tot}` or the Ohmic current density 
:math:`j_\Omega` that is known. By default, it is assumed that it is the total 
that is provided, from which the Ohmic component is then calculated by subtracting
the bootstrap current density (and other components e.g. the runaway current).

However, if it is the Ohmic component that is known, simply include the following
in the simulation setup script:

.. code-block:: python

   ds.eqsys.j_bs.setInitMode(BootstrapCurrent.BOOTSTRAP_INIT_MODE_OHMIC)


Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.BootstrapCurrent.BootstrapCurrent
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__
