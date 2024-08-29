.. _ds-eqsys-E_field:

ElectricField
=============
The ``ElectricField`` class holds settings for the parallel electric field solved for by DREAM, 
which is described by the unknown ``E_field``. This quantity is defined as

.. math::
   \mathrm{E\_field} = \frac{\langle \boldsymbol{E}\cdot\boldsymbol{B} \rangle}{\sqrt{\langle B^2 \rangle}},

where angle brackets :math:`\langle\cdot\rangle` denote a flux-surface average, 
:math:`\boldsymbol{E}` denotes the electric field, :math:`\boldsymbol{B}` the 
magnetic field and :math:`\varphi` the toroidal angle.

The electric field can be solved for
in three different ways:

(1) By prescribing the electric field profile in time :math:`E_\parallel = \tilde{E}(t,r)` (``TYPE_PRESCRIBED``)
(2) By solving the induction equation (``TYPE_SELFCONSISTENT``):
(3) By requiring Ohm's law to be satisfied, :math:`j_\Omega=\sigma E_\parallel`,
    and providing the ohmic current density profile :math:`j_\Omega`, along with
    the temperature and density required to determine the plasma conductivity
    :math:`\sigma`.

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

   where :math:`\psi_{\rm p}` is the poloidal flux, :math:`\mu_0` 
   is the vacuum permeability, :math:`V'` the
   spatial Jacobian and :math:`r` the minor radius in the outer midplane. 
   The loop voltage :math:`V_{\rm loop}` is defined by

   .. math::
      V_{\rm loop} = 2\pi \frac{\langle \boldsymbol{E}\cdot\boldsymbol{B} \rangle}{\langle \nabla \varphi \cdot \boldsymbol{B} \rangle}.

.. contents:: Page overview
   :local:
   :depth: 3

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

   ds.eqsys.E_field.setInitialProfile(efield=E_field, radius=E_field_r)

where ``E_field`` is a vector of size ``nr`` giving the initial electric field
profile, and ``E_field_r`` is a vector representing the radial grid on which the
initial electrc field profile is defined.

Alternatively, the electric field can be initialized such that it gives the
desired ohmic current density profile (achieved by inverting Ohm's law):

.. code-block:: python

   ds.eqsys.j_ohm.setInitialProfile(j=j_ohm, radius=j_ohm_r, Ip0=Ip0)

where ``j_ohm`` is the desired current density profile and ``j_ohm_r`` its grid.
If the optional parameter ``Ip0`` is provided, the given current density
profile is rescaled such that the prescribed ohmic current density profile
yields exactly a total plasma current ``Ip0`` when integrated on the internal
DREAM grid (which is not necessarily the same as the grid on which ``j_ohm`` is
given).

Boundary condition
******************
The second equation in :eq:`eq_E_field` requires a boundary condition at
:math:`r=r_{\rm wall}` to be given. In DREAM, three different boundary
conditions can be applied at the tokamak wall.

**The first boundary condition**, ``BC_TYPE_PRESCRIBED``, prescribes the time
evolution of the loop voltage :math:`V_{\rm loop,wall}` on the tokamak wall
(normalized to the tokamak major radius :math:`R_0`). This is particularly
useful for simulating experimental scenarios where the parameter
:math:`V_{\rm loop,wall}` has been measured.

.. code-block:: python

   import DREAM.Settings.Equations.ElectricField as ElectricField

   ds = DREAMSettings()
   ...
   # Tokamak major radius
   R0        = 1.65
   # Define evolution of V_loop/R0
   Vmax      = 1
   V_loop_t  = np.linspace(0, 1, 100)
   V_loop_R0 = (Vmax/R0)*(1 - (1-t)**2)

   ds.eqsys.E_field.setBoundaryCondition(bctype=ElectricField.BC_TYPE_PRESCRIBED,
                                         V_loop_wall_R0=V_loop_R0, times=V_loop_t, R0=R0)

**The second boundary condition**, ``BC_TYPE_SELFCONSISTENT``, instead lets the user
specify the (inverse) tokamak wall time, :math:`1/\tau_{\rm wall}`, which
directly corresponds to the wall resistivity. The perfectly conducting limit
:math:`\tau_{\rm wall} = \infty` is supported, and is obtained by setting
``inverse_wall_time = 0``. In case a cylindrical geometry is used, the major
radius ``R0`` can be explicitly set, independently of the geometry used, since 
the external inductance otherwise diverges for infinite major radius.

.. code-block:: python

   import DREAM.Settings.Equations.ElectricField as ElectricField

   ds = DREAMSettings()
   ...
   # Tokamak major radius
   R0       = 1.65
   # Wall time
   tau_wall = .01   # (s)
   
   ds.eqsys.E_field.setBoundaryCondition(bctype=ElectricField.BC_TYPE_SELFCONSISTENT,
                                         inverse_wall_time=1/tau_wall, R0=R0)

**The third boundary condition**, ``BC_TYPE_TRANSFORMER``, can be seen as a
combination of the two other boundary conditions, with a resistive wall *and*
a prescribed loop voltage (although at the transformer, which passes through the
axis of symmetry at :math:`R=0`).

.. code-block:: python

   import DREAM.Settings.Equations.ElectricField as ElectricField

   ds = DREAMSettings()
   ...
   # Tokamak major radius
   R0        = 1.65
   # Define evolution of V_loop/R0
   Vmax      = 1
   V_loop_t  = np.linspace(0, 1, 100)
   V_loop_R0 = (Vmax/R0)*(1 - (1-t)**2)
   # Wall time
   tau_wal   = .01  # (s)

   ds.eqsys.E_field.setBoundaryCondition(bctype=ElectricField.BC_TYPE_PRESCRIBED,
                                         inverse_wall_time=1/tau_wall, R0=R0,
                                         V_loop_wall_R0=V_loop_R0, times=V_loop_t)

.. note::

   With all boundary conditions, the tokamak wall location ``wall_radius`` must
   be specified. This parameter denotes the minor radial coordinate of the
   tokamak wall (i.e. the distance of the wall from the center of the plasma).
   This is done via the call ``ds.radialgrid.setWallRadius(b)``, where ``b`` is
   the desired wall radius.


Prescribed ohmic current
------------------------
If the desired ohmic current density profile :math:`j_\Omega` is known, the
corresponding electric field can be found by inverting Ohm's law

.. math::

   \frac{j_\Omega}{B} = \sigma \frac{\left\langle\boldsymbol{E}\cdot\boldsymbol{B}\right\rangle}{\left\langle B^2 \right\rangle}

where :math:`\sigma` denotes the plasma conductivity, which depends on primarily
the electron density and temperature. To evaluate the electric field from this
law, the following configuration can be made:

.. code-block:: python

   import DREAM.Settings.Equations.ElectricField as ElectricField

   ds.eqsys.E_field.setType(ElectricField.TYPE_PRESCRIBED_OHMIC_CURRENT)
   ds.eqsys.j_ohm.setCurrentProfile(j=j_ohm, radius=r, times=t, Ip0=Ip0)

where ``j_ohm`` is the desired current density profile, ``r`` is the radial
grid on which the profile is specified, and ``t`` is its time grid. If no time
evolving current density is needed, ``j_ohm`` can be given as a 1D array. If the
optional scalar parameter ``Ip0`` is given, the current density profile is
rescaled such that, when integrated over the plasma cross-section, it yields
exactly the total plasma current ``Ip0``.

Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.ElectricField.ElectricField
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __getitem__

