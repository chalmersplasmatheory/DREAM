.. _ds-eqsys-E_field:

ElectricField
=============
The ``ElectricField`` class holds settings for the parallel electric field
:math:`E_\parallel` solved for by DREAM. The electric field can be solved for
in two different ways:

(1) By prescribing the electric field profile in time :math:`E_\parallel = \tilde{E}(t,r)` (``TYPE_PRESCRIBED``)
(2) By solving the induction equation (``TYPE_SELFCONSISTENT``):

   .. math::
      :label: eq_E_field

       \begin{cases}
           \frac{\partial\psi_{\rm p}}{\partial t} &= V_{\rm loop},\\
           \mu_0\frac{j_\parallel}{B}\left\langle \boldsymbol{B}\cdot\nabla\varphi \right\rangle &=
           \frac{1}{V'}\frac{\partial}{\partial r}\left[
               V'\left\langle \frac{\left|\nabla r\right|^2}{R^2} \right\rangle
               \frac{\partial\psi_{\rm p}}{\partial r}
           \right]
       \end{cases}

   where :math:`\psi_{\rm p}` is the poloidal flux, :math:`V_{\rm loop}\equiv 2\pi RE_\parallel`
   the loop voltage (at major radius :math:`R`), :math:`\mu_0` is the vacuum permeability,
   :math:`\boldsymbol{B}` the magnetic field, :math:`\varphi` the toroidal angle, :math:`V'` the
   spatial Jacobian, and :math:`r` the minor radius in the outer midplane. Angle brackets
   :math:`\langle\cdot\rangle` denote a flux-surface average.

Prescribed evolution
--------------------
The electric field can be prescribed in both space and time. The full syntax for
prescribing an electric field evolution is

.. code-block:: python

   ds.eqsys.E_field.setPrescribedData(efield=E_field, radius=E_field_r, times=E_field_t)

where ``E_field`` is an array of shape ``(nt, nr)`` specifying the parallel
electric field strength in space and time, ``E_field_r`` is a vector with ``nr``
elements representing the radial grid on which the electric field is prescribed,
and ``E_field_t`` is a vector with ``nt`` elements representing the time grid on
which the electric field is prescribed.

For convenience, it is possible to prescribe a constant and radially uniform 
electric field by only specifying it as a scalar:

.. code-block:: python

   ds.eqsys.E_field.setPrescribedData(0.3)

This will cause the electric field to be :math:`E_\parallel = 0.3\,\mathrm{V/m}`
at every radius, in every time step, during the simulation.

Self-consistent evolution
-------------------------
In this mode the electric field is evolved by solving the system
:eq:`eq_E_field`. The evolution of the electric field is thus coupled to the
poloidal flux :math:`\psi_{\rm p}` and, by extension, to the total plasma
current density :math:`j_{\rm tot}`. To evolve the electric field
self-consistently, the first thing one must do is to set the type of the
electric field to ``TYPE_SELFCONSISTENT`` via a call to
:py:meth:`~DREAM.Settings.Equations.ElectricField.setType`:

.. code-block:: python

   import DREAM.Settings.Equations.ElectricField as ElectricField

   ...

   ds.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)

By default, the initial electric field profile will be identically zero. This
is automatically overridden if the settings are loaded from an old output file
(see :ref:`restart`), but the initial profile can also be set explicitly via a
call to :py:meth:`~DREAM.Settings.Equations.ElectricField.setInitialProfile`:

.. code-block:: python

   ds.eqsys.setInitialProfile(efield=E_field, radius=E_field_r)

where ``E_field`` is a vector of size ``nr`` giving the initial electric field
profile, and ``E_field_r`` is a vector representing the radial grid on which the
initial electrc field profile is defined.

Boundary condition
******************
The second equation in :eq:`eq_E_field` requires a boundary condition at
:math:`r=r_{\rm wall}` to be given. In DREAM, two different boundary conditions
can be applied at the tokamak wall.

The first boundary condition, ``BC_TYPE_PRESCRIBED``, prescribes the time
evolution of the loop voltage :math:`V_{\rm loop,wall}` on the tokamak wall.
This is particularly useful for simulating experimental scenarios where the
parameter :math:`V_{\rm loop,wall}` has been measured.

The second boundary condition, ``BC_TYPE_SELFCONSISTENT``, instead lets the user
specify the (inverse) tokamak wall time, :math:`1/\tau_{\rm wall}`, which
directly corresponds to the wall resistivity. The perfectly conducting limit
:math:`\tau_{\rm wall} = \infty` is supported, and is obtained by setting
``inverse_wall_time = 0``.

.. note::

   With both boundary conditions, the tokamak wall location ``wall_radius`` must
   be specified. This parameter denotes the minor radial coordinate of the
   tokamak wall (i.e. the distance of the wall from the center of the plasma).

Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.ElectricField.ElectricField
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __getitem__


