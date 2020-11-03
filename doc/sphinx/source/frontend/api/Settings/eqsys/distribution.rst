.. _ds-eqsys-distfunc:

DistributionFunction
====================
The ``DistributionFunction`` class contains all settings for the hot and runaway
electron distribution functions. The distribution functions are automatically
enabled whenever their corresponding grids (:ref:`hottailgrid<ds-momentumgrid>`
and :ref:`runawaygrid<ds-momentumgrid>`) are enabled.

The distribution functions are solved for using the 3D Fokker-Planck equation

.. math::

   \frac{\partial f}{\partial t} +
   eE\left[
       \frac{1}{p^2}\frac{\partial}{\partial p}\left( p^2\xi f \right) +
       \frac{1}{p}\frac{\partial}{\partial\xi}\left[ \left( 1-\xi^2 \right) f \right]
   \right] =
   C\left\{ f \right\} +
   \frac{1}{\mathcal{V}'}\frac{\partial}{\partial r} \left[
       \mathcal{V}'\left( -Af + D\frac{\partial f}{\partial r} \right)
   \right] + S,

where the second term on the LHS represents electric field acceleration along
magnetic field lines, :math:`C\{f\}` denotes the Fokker-Planck collision
operator, the second term on the RHS represents advective-diffusive radial
transport with arbitrary advection and diffusion coefficients :math:`A` and
:math:`D`, and :math:`S` is a collection of source terms (most notably the
production of secondary runaway electrons through the avalanche mechanism
and the particle source term which ensures quasi-charge neutrality).
Additional terms not shown here can also be added, as described further down on
this page.

For the hot electron distribution function, it is also possible to solve the
pitch angle-averaged 2D Fokker-Planck equation by setting ``nxi=1`` on the
:ref:`hottailgrid<ds-momentumgrid>`. In this case, the equation is reduced
with the ordering :math:`\nu_D \gg eE/p \gg \mathrm{remainder}`, in which case the 
pitch-angle distribution takes an analytic form with which the equation is 
analytically averaged. The approximation becomes increasingly better for higher 
plasma effective charge and lower energy. It is inadequate for describing
Dreicer seed generation or the shape of the runaway tail, but can provide a 
fair approximation for the runaway rate of other mechanisms. 


.. contents:: Page overview
   :local:
   :depth: 3

.. _advection-interpolation:

Advection interpolation
***********************
In the finite volume method, which is used for discretizing most equations in
DREAM (particularly the Fokker-Planck equation), advection terms are discretized
as

.. math::

   \left[ \frac{1}{\mathcal{V}'}\frac{\partial}{\partial z^k}\left( \mathcal{V}' A_k f \right) \right]_i =
   \frac{\mathcal{V}_{i+1/2}' A_{k,i+1/2} f_{i+1/2} - \mathcal{V}_{i-1/2}' A_{k,i-1/2} f_{i-1/2}}
   {\mathcal{V}_i' \Delta z^k}

where half-indices are located on cell faces, while the distribution function
:math:`f` is given at cell centres, corresponding to integer indices. Thus, to
obtain a practical discretization for the advection term it becomes necessary to
interpolate in :math:`f`. This is generally achieved by introducing a set of
interpolation coefficients :math:`\left\{ \delta^{(i)}_j \right\}` and letting

.. math::

   f_{i-1/2} = \sum_{j=-2}^1 \delta^{(i)}_j f_{i+j}.

The interpolation coefficients can be chosen in a number of ways, and one of the
simplest choices imaginable is :math:`\delta^{(i)}_{-2} = \delta^{(i)}_{1} = 0`,
:math:`\delta^{(i)}_{-1} = \delta^{(i)}_{0} = 1/2`. This is the so-called
``CENTRED`` interpolation scheme, which is used by default in DREAM. For
strongly advective problems, however, it is usually desirable to choose these
coefficients more carefully.

Advective problems can famously yield unstable solutions containing spurious,
unphysical oscillations. These oscillations result from the unphysical nature of
the interpolation done with schemes such as the ``CENTRED`` described above and
can in fact be avoided with other interpolation schemes. DREAM provides a number
of such *positivity-preserving* interpolation schemes, which avoid oscillations
and ensure that the resulting distribution function is everywhere positive. The
simplest schemes are the linear first and second-order upwind schemes, which
are computationally cheap but less accurate. For most cases, the non-linear
*flux limiter* schemes are to be preferred.

.. note::

   All flux limiter schemes require the :ref:`non-linear solver<ds-solver>` to
   be used.

List of schemes
---------------

The following interpolation schemes are provided by DREAM:

+--------------------------------+-------------+------------------------------------------------------------------------------------+
| Name                           | Non-linear? | Reference                                                                          |
+================================+=============+====================================================================================+
| ``AD_INTERP_CENTRED``          | No          |                                                                                    |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_UPWIND``           | No          |                                                                                    |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_UPWIND_2ND_ORDER`` | No          |                                                                                    |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_DOWNWIND``         | No          |                                                                                    |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_QUICK``            | No          | `Leonard (1980) <https://ui.adsabs.harvard.edu/abs/1980cmf..book..159L/abstract>`_ |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_MUSCL``            | **Yes**     | `van Leer (1979) <https://doi.org/10.1016/0021-9991(79)90145-1>`_                  |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_SMART``            | **Yes**     | `Gaskell and Lau (1988) <https://doi.org/10.1002/fld.1650080602>`_                 |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_OSPRE``            | **Yes**     | Waterson and Deconinck (1995)                                                      |
+--------------------------------+-------------+------------------------------------------------------------------------------------+
| ``AD_INTERP_TCDF``             | **Yes**     | `Zhang et al (2015) <https://doi.org/10.1016/j.jcp.2015.08.042>`_                  |
+--------------------------------+-------------+------------------------------------------------------------------------------------+

Two useful systematic reviews of advection interpolation schemes and flux limiter methods
are given by `N P Waterson & H Deconinck (2007) <https://doi.org/10.1016/j.jcp.2007.01.021>`_ 
and `Zhang et al (2015) <https://doi.org/10.1016/j.jcp.2015.08.042>`_.

.. note::

   For most applications where flux limiters are desired, we recommend the TCDF
   flux limiter scheme.

For flux limiter schemes, it is also possible to specify a *flux limiter
damping* parameter. This parameter, here denoted :math:`\eta`, is used when
updating the interpolation coefficients in the non-linear iteration process to
under-relax the algorithm as this may accelerate convergence. In each iteration
the new interpolation coefficient is set to

.. math::

   \delta = \delta_{\rm prev} + \eta\left( \delta_{\rm calc} - \delta_{\rm prev} \right),

where :math:`\delta_{\rm prev}` is the old value of the interpolation
coefficient and :math:`\delta_{\rm calc}` is the calculated new value.

Advection jacobian mode
^^^^^^^^^^^^^^^^^^^^^^^
The flux limiter interpolation schemes are non-linear as the coefficients depend 
on the gradients of the unknown in which they interpolate. There are three different
methods of treating the jacobian of the interpolation coefficients:

+-------------------------------+--------------------------------------------------------------------------+
| Name                          | Description                                                              |
+===============================+==========================================================================+
| ``AD_INTERP_JACOBIAN_LINEAR`` | Default. Neglects any contribution from the unknown to the jacobian      |
+-------------------------------+--------------------------------------------------------------------------+
| ``AD_INTERP_JACOBIAN_FULL``   | Fully accounts for the unknown contribution to the jacobian              |
+-------------------------------+--------------------------------------------------------------------------+
| ``AD_INTERP_JACOBIAN_UPWIND`` | Uses the jacobian corresponding to ``UPWIND`` interpolation coefficients |
+-------------------------------+--------------------------------------------------------------------------+

.. note:: 
   ``FULL`` is only active for the smooth ``TCDF`` and ``OSPRE`` flux limiter 
   schemes, for all others it defaults to ``LINEAR``. 
   ``UPWIND`` requires a non-linear solver as the jacobian does not match the
   function vector, and it will sometimes take many iterations to converge,
   but each iteration is faster due to fewer non-zero elements in the matrix. 


Example
-------
The following example illustrates how to use the TCDF flux limiter with full 
jacobian in all dimensions on the hot electron grid:

.. code-block:: python

   import DREAM.Settings.Equations.DistributionFunction as DistFunc

   ds = DREAMSettings()
   ...
   ds.eqsys.f_hot.setAdvectionInterpolationMethod(
      ad_int=DistFunc.AD_INTERP_TCDF,
      ad_jac=DistFunc.AD_INTERP_JACOBIAN_FULL)

It is also possible to specify different interpolation methods in different
dimensions. Using a flux limiter scheme may be more important in the :math:`p`
dimension due to the strong advective flows due to the electric field, whereas
the pitch dimension may be dominated by (diffusive) pitch-angle scattering, and
the radial dimension has no fluxes at all. In this case, it may be sufficient to
apply the flux limiter to the :math:`p` dimension only:

.. code-block:: python

   import DREAM.Settings.Equations.DistributionFunction as DistFunc

   ds = DREAMSettings()
   ...
   ds.eqsys.f_hot.setAdvectionInterpolationMethod(ad_int_p1=DistFunc.AD_INTERP_TCDF)

The other dimensions default to the ``AD_INTERP_CENTRED`` scheme.

.. warning::

   Although flux limiters will generally produce smooth and stable solutions,
   they are **NOT** at all guaranteed to be accurate. You should always verify
   that the results you obtain are correct by varying the various resolution
   parameters of your simulation (usually ``nr``, ``np``, ``nxi`` and ``nt``).
   Convergence of the non-linear solver is not guaranteed with flux limiters,
   as sometimes there may be two valid solutions that the system oscillates 
   between. It can help to increase resolution, try ``fluxlimiterdamping``,
   another jacobian mode or another flux limiter. 

Boundary conditions at pMax
***************************
Numerical solution of the Fokker-Planck equation on a finite momentum grid
requires boundary condition to specified at :math:`p=p_{\rm max}`, the upper
boundary in the momentum dimension. DREAM supports three different forms for
this boundary condition, corresponding to

(1) setting :math:`f(p>p_{\rm max}) = 0` (default)
(2) ensuring that the flux out across :math:`p=p_{\rm max}` is the same as the
    flux into the last cell on the momentum grid.
(3) extrapolating the flux out across :math:`p=p_{\rm max}` based on the flux
    into the two last cells on the momentum grid.

**Which boundary condition should I use?**
You can choose whichever of the boundary conditions you prefer. The conditions
(1) and (2) may be marginally cheaper to use, but all are equally arbitrary. If
the boundary condition appears to affect the result of your simulations, you
should rather increase :math:`p_{\rm max}` so that it is located sufficiently
far from where the interesting physics are happening.

In some cases, the boundary condition may give rise to numerical instabilities
on the boundary. In these cases it is well motivated to try an alternative
boundary condition. If none of the boundary conditions are able to stabilize
the solution near the boundary, you should proceed to applying a
:ref:`flux limiter scheme<advection-interpolation>`.

List of boundary conditions
---------------------------

+-------------------+------------------------------------------------------------------------------------------------------------------+
| Name              | Description                                                                                                      |
+===================+==================================================================================================================+
| ``BC_F_0``        | Assume :math:`f(p>p_{\rm max}) = 0` (default)                                                                    |
+-------------------+------------------------------------------------------------------------------------------------------------------+
| ``BC_PHI_CONST``  | Assume :math:`\Phi^{(p)}_{N_p+1/2} = \Phi^{(p)}_{N_p-1/2}`                                                       |
+-------------------+------------------------------------------------------------------------------------------------------------------+
| ``BC_DPHI_CONST`` | Assume :math:`\Phi^{(p)}_{N_p+1/2} = \left( 1 - \delta \right)\Phi^{(p)}_{N_p-1/2} + \delta\Phi^{(p)}_{N_p-3/2}` |
+-------------------+------------------------------------------------------------------------------------------------------------------+

.. note::

   When using both hot and runaway electron grids, the boundary condition at
   :math:`p_{\rm hot} = p_{\rm hot,max}` is automatically set to allow particles
   to flow freely across the boundary between the grids. It is then only
   necessary to specify the boundary condition at the
   :math:`p_{\rm RE} = p_{\rm RE,max}` boundary on the runaway grid.

.. note::

   All particles leaving the hot electron grid through :math:`p=p_{\rm max}`
   are automatically moved to the runaway density :ref:`n_re<ds-eqsys-n_re>`.
   When particles leave the runaway electron grid through :math:`p=p_{\rm max}`
   they are still counted as part of :ref:`n_re<ds-eqsys-n_re>`, and thus
   this quantity may deviate from the density moment of ``f_re``.

Example
-------
The following example shows how to use the ``BC_PHI_CONST``:

.. code-block:: python

   import DREAM.Settings.Equations.DistributionFunction as DistFunc

   ds = DREAMSettings()
   ...
   ds.eqsys.f_hot.setBoundaryCondition(DistFunc.BC_PHI_CONST)


Initialization
**************
An initial condition is required when solving the Fokker-Planck equation and in
DREAM this can be either in the form of a Maxwellian with prescribed density and
temperature profiles, or as an arbitrary prescribed function.

As Maxwellian
-------------
DREAM can initialize the distribution function to a Maxwell-Jüttner distribution
function:

.. math::

   f(r,p,\xi) = \frac{n_0(r)}{4\pi\Theta_0(r) K_2\left(1/\Theta_0(r)\right)} e^{-\gamma/\Theta_0(r)},

where :math:`n_0(r)` is the initial electron density profile,
:math:`\Theta_0(r) = T_0(r)/m_ec^2` is the initial electron temperature profile
normalized to the electron rest mass, :math:`\gamma` is the Lorentz factor and
:math:`K_2(x)` is the second-order modified Bessel function of the second kind.

To prescribe the initial distribution function as a Maxwellian, the user should
call the method ``setInitialProfiles()``, providing the density and temperature
profiles determining the shape of the initial distribution function at every
radius.

Example
^^^^^^^
The following example illustrates how to prescribe the initial hot electron
distribution function as a Maxwellian:

.. code-block:: python

   ds = DREAMSettings()
   ...
   # Radial grid (may be different for n0 and T0)
   r = np.linspace(r0, r1, nr)

   # Density profile
   n0 = np.array([...])     # shape (nr,)  [m^-3]
   # Temperature profile
   T0 = np.array([...])     # shape (nr,)  [eV]

   ds.eqsys.f_hot.setInitialProfiles(n0=n0, T0=T0, rn0=r, rT0=r)

If uniform density and/or temperature profiles are desired, ``n0`` and/or
``T0`` can be given as scalars and their corresponding radial grid parameters
omitted from the call to ``setInitialProfiles()``:

.. code-block:: python

   ds = DREAMSettings()
   ...
   # Uniform density of n=2x10^19 m^-3 and temperature T=1 keV
   ds.eqsys.f_hot.setInitialProfiles(n0=2e19, T0=1000)


Arbitrary distribution function
-------------------------------
Both ``f_hot`` and ``f_re`` can also be initialized with arbitrary distribution
functions. The functions should be represented as 3D arrays with the shape
``(nr, np2, np1)`` (for :math:`p/\xi` coordinates, ``np1=np`` and ``np2=nxi``).
Initialization is done via a call to the method ``setInitialValue()``.

Example
^^^^^^^
The following example illustrates how to prescribe an arbitrary initial hot
electron distribution function:

.. code-block:: python

   ds = DREAMSettings()
   ...
   r  = np.linspace(r0, r1, nr)
   p  = np.linspace(p0, p1, np)
   xi = np.linspace(xi0, xi1, nxi)

   f0 = np.array([...])     # shape (nr, nxi, np)

   ds.eqsys.f_hot.setInitialValue(f=f0, r=r, p=p, xi=xi)


Magnetic ripple
***************
Due to the finite number of toroidal field coils used for tokamaks, the toroidal
symmetry of the devices is somewhat broken. The perturbed magnetic field can
lead to both increased particle transport, as well as increased pitch-angle
scattering.

In DREAM, we can account for the increased pitch-angle scattering using the
diffusion operator derived by
`Kurzan et al (1995) <https://doi.org/10.1103/PhysRevLett.75.4626>`_, which for DREAM
takes the form

.. math::

   D^{\xi_0\xi_0} =
   \frac{\pi}{32}\frac{eB}{m_e}
   \frac{B_{\rm min}}{B}\frac{\xi^2}{\xi_0^2}\frac{1-\xi_0^2}{p^2}
   \frac{v_\parallel}{c}\left( \frac{\delta B_{nm}}{B} \right)^2
   H_{nm}(p\xi),

with

.. math::

   H_{nm}(p\xi) = \begin{cases}
       \frac{1}{\Delta p_{nm}}, \qquad & \left|p_\parallel - p_{nm} \right| < \Delta p_{nm},\\
       0, \qquad & \text{otherwise}
   \end{cases}

and resonant momentum and resonance width

.. math::

   p_{nm} &= \frac{e G(r)}{m_e cn N_{\rm c}},\\
   \Delta p_{nm} &= \sqrt{\frac{\delta B_{nm}}{B} p_\parallel p_\perp}.

where :math:`G(r)` gives the toroidal magnetic field variation and
:math:`N_{\rm c}` denotes the number of toroidal field coils. Since the
derivation assumes small pitch angles and slab geometry, the geometric factor
:math:`\xi^2B_{\rm min}/\xi_0^2B = 1`.

Input data
----------
To include the magnetic ripple pitch scattering in simulation, it is necessary
to provide the number of magnetic field coils :math:`N_{\rm c}`, as well as the
perturbation magnitude :math:`\delta B_{nm}/B` at one or more mode numbers
:math:`(n,m)` (toroidal, poloidal). The magnetic perturbations can be specified
as functions of time and radius, and will be interpolated as necessary.

In cylindrical (radial grid) mode, the function :math:`G(r)` becomes infinite.
However, this is compensated by the fact that a large-aspect ratio tokamak would
have an infinite number of toroidal field coils, :math:`N_{\rm c}`. To still
simulate the ripple in a cylindrical field, we therefore introduce the distance
between magnetic field coils

.. math::

   \Delta_{\rm coils} = \frac{2\pi R_0}{N_{\rm c}}.

The resonant momentum then takes the form

.. math::

   p_{nm} = \frac{e G(r) \Delta_{\rm coils}}{2\pi m_ec n R_0},

with :math:`G(r)/R_0` forming the finite toroidal magnetic field stregth.

Resonant momentum
-----------------
When including the ripple pitch scattering in a simulation, it is possible to
save the resonance momentum :math:`p_{nm}` calculated at each time step, and for
each pair of mode numbers, as an :ref:`other quantity<otherquantity>`. The
available ripple quantities are

+----------------------------+--------------------------------------------------------+
| Name                       | Description                                            |
+============================+========================================================+
| ``fluid/f_hot_ripple_m``   | Poloidal mode numbers corresponding to each resonance. |
+----------------------------+--------------------------------------------------------+
| ``fluid/f_hot_ripple_n``   | Toroidal mode numbers corresponding to each resonance. |
+----------------------------+--------------------------------------------------------+
| ``fluid/f_hot_ripple_pmn`` | Resonant momentum for ``f_hot``.                       |
+----------------------------+--------------------------------------------------------+
| ``fluid/f_re_ripple_m``    | Poloidal mode numbers corresponding to each resonance. |
+----------------------------+--------------------------------------------------------+
| ``fluid/f_re_ripple_n``    | Toroidal mode numbers corresponding to each resonance. |
+----------------------------+--------------------------------------------------------+
| ``fluid/f_re_ripple_pmn``  | Resonant momentum for ``f_re``.                        |
+----------------------------+--------------------------------------------------------+

Example
-------
This example shows how to include the magnetic ripple pitch scattering mechanism
in a simulation of the hot electron distribution function:

.. code-block:: python

   ds = DREAMSettings()
   ...
   # Number of toroidal field coils
   nCoils = 12

   # Include three resonances
   dB_B   = [1e-4,5e-5,2e-5]
   m      = [1,1,1]
   n      = [1,2,3]

   # Apply settings
   ds.eqsys.f_hot.setRipple(m=m, n=n, dB_B=dB_B, ncoils=nCoils)

In cylindrical grid mode, the parameter ``deltaCoil`` should be specified
instead of ``ncoils``.

The arrays ``m`` and ``n`` should always have the same number of elements, which
should also be the same number of elements as in the first dimension of ``dB_B``.

To prescribe a set of spatiotemporally varying perturbations, replace the
perturbations in the example above with 2D arrays and specify time and radial
grid vectors:

.. code-block:: python

   import numpy as np

   ds = DREAMSettings()

   t = np.linspace(0, 1, 10)
   r = np.linspace(0, 0.5, 20)

   dB_B_11 = np.array([...])
   dB_B_21 = np.array([...])
   dB_B_31 = np.array([...])

   dB_B = [dB_B_11, dB_B_21, dB_B_31]
   m    = [1,1,1]
   n    = [1,2,3]

   ds.eqsys.f_hot.setRipple(m=m, n=n, dB_B=dB_B, r=r, t=t, ncoils=nCoils)


Radial Transport
****************
A radial transport term can be added to the Fokker-Planck equation. The general
transport term takes the form

.. math::

   \nabla\cdot\boldsymbol{\Phi}_{\rm transport} = \frac{1}{\mathcal{V}'}\frac{\partial}{\partial r}\left[
       \mathcal{V}'\left( A f + D\frac{\partial f}{\partial r} \right)
   \right],

where :math:`A` and :math:`D` are the advection and diffusion coefficients
respectively. DREAM provides two different ways of prescribing these
coefficients:

(1) Let the user prescribe :math:`A` and/or :math:`D` directly.
(2) Use the model derived by Rechester and Rosenbluth and let the user prescribe
    the magnetic perturbation level :math:`\delta B/B`.


Prescribed coefficients
-----------------------
The user can directly prescribe the coefficients :math:`A=A(t,r,p,\xi)` and
:math:`D=D(t,r,p,\xi)` by calling the methods ``transport.prescribeAdvection()``
and ``transport.prescribeDiffusion()``. The coefficients can be either scalars
(in which case they are assumed to be constant and uniform in all parameters) or
4D arrays, representing the coefficient's dependence on time, radius, momentum
and pitch, in that order.

Example
^^^^^^^
The following example shows how to prescribe advective-diffusive transport for
the hot electron distribution function:

.. code-block:: python

   ds = DREAMSettings()
   ...
   t  = np.linspace(t0, t1, nt)
   r  = np.linspace(r0, r1, nr)
   p  = np.linspace(p0, p1, np)
   xi = np.linspace(xi0, xi1, nxi)

   A = np.array([...])  # shape (nt, nr, np, nxi)
   D = np.array([...])  # shape (nt, nr, np, nxi)

   ds.eqsys.f_hot.transport.prescribeAdvection(ar=A, t=t, r=r, p=p, xi=xi)
   ds.eqsys.f_hot.transport.prescribeDiffusion(drr=D, t=t, r=r, p=p, xi=xi)

Contrary to what the example may seem to suggest, the coefficients do *not* need
to share the same time, radius, momentum and pitch grids. The grids are also
*not* required to be the same as those used internally by DREAM for the
simulation, as the provided coefficients will be interpolated as needed.

If a constant coefficient is desired, the following more compact syntax may be
used:

.. code-block:: python

   ds.eqsys.f_hot.transport.prescribeDiffusion(10)

which will prescribe a constant diffusion coefficient with the value
:math:`10\,\mathrm{m}^2/\mathrm{s}`.

Rechester-Rosenbluth
--------------------
Rechester and Rosenbluth derived a diffusion operator for describing the radial
transport due to turbulent perturbations of the magnetic field. Given a magnetic
field perturbation strength :math:`\delta B/B`, the diffusion coefficient for
the Rechester-Rosenbluth operator takes the form

.. math::

   D = \pi q R_0\left(\frac{\delta B}{B}\right)^2 \left|v_\parallel\right|,

where :math:`q(r)` is the tokamak safety factor, :math:`R_0` is the tokamak
major radius and :math:`v_\parallel` is the particle speed along the magnetic
field line. In DREAM, the user provides the magnetic perturbation
:math:`\delta B/B` as a function of time and radius.

.. note::

   At the time of writing, DREAM does not have a simple way of determining the
   safety factor :math:`q(r)` for any radial grid. Because of this, the safety
   factor is set to :math:`q(r)\equiv 1` in the implementation of this operator,
   effectively requiring the user to provide this dependence in the input
   parameter :math:`\delta B/B`.

Example
^^^^^^^
The following example illustrates how to prescribe a Rechester-Rosenbluth
diffusion coefficient in DREAM:

.. code-block:: python

   ds = DREAMSettings()
   ...
   # Grids on which the perturbation is defined
   t = np.linspace(t0, t1, nt)
   r = np.linspace(r0, r1, nr)

   # Magnetic perturbation
   dB_B = np.array([...])   # shape (nt, nr)

   ds.eqsys.f_hot.setMagneticPerturbation(dBB=dB_B, t=t, r=r)

If a constant and uniform perturbation level is desired, the more compact syntax

.. code-block:: python

   ds.eqsys.f_hot.setMagneticPerturbation(dBB=1e-3)

can also be used. In this example, the magnetic perturbation
:math:`\delta B/B = 10^{-3}` is used in every time step and at all radial
positions.


Boundary conditions
-------------------
The general transport term requires a boundary condition to specified at the
plasma boundary :math:`r=a`. DREAM currently supports two different boundary
conditions at this boundary:

(1) Assume that no particles can cross the plasma boundary (default).
(2) Assume that :math:`f(r>a)=0` and that fluxes are set accordingly.

Condition (1) is a reflective boundary condition which naturally conserves the
total number of particles in :math:`f`. It prevents particles from leaving the
plasma.

Condition (2) generally results in a net outflux of particles, and allows
particles to be lost from the system. However, to conserve quasi-neutrality,
DREAM compensates for the lost electrons by adding additional cold electrons.

List of options
^^^^^^^^^^^^^^^
The following boundary conditions can be applied:

+---------------------+------------------------------------------------------------------------------+
| Name                | Description                                                                  |
+=====================+==============================================================================+
| ``BC_CONSERVATIVE`` | Assume that no particles can leave the plasma through :math:`r=a`. (default) |
+---------------------+------------------------------------------------------------------------------+
| ``BC_F_0``          | Assume that :math:`f(r>a) = 0`.                                              |
+---------------------+------------------------------------------------------------------------------+

Example
^^^^^^^
The following example illustrates how to apply a lossy boundary condition to
the hot electron distribution function:

.. code-block:: python

   import DREAM.Settings.TransportSettings as Transport

   ds = DREAMSettings()
   ...
   ds.eqsys.f_hot.transport.setBoundaryCondition(Transport.BC_F_0)


Runaway avalanching
*******************
Secondary runaway electrons (produced through the so-called avalanche mechanism)
can be included in the simulation of the distribution function through a
Rosenbluth-Putvinski source term
(`Rosenbluth and Putvinski 1997 <https://doi.org/10.1088/0029-5515/37/10/I03>`_).
This source term is however enabled on the :ref:`n_re <ds-eqsys-n_re>` settings
object using the ``setAvalanche()`` method with the option
``RunawayElectrons.AVALANCHE_MODE_KINETIC``.

Since the source term diverges in the origin, creating an infinite number of knock-ons,
a cutoff parameter ``pCutAvalanche`` is introduced, defined as the minimum momentum
to which the source is applied. 

.. warning::
   The momentum cutoff must be chosen to be smaller than the critical momentum for 
   runaway generation, otherwise the cutoff will act as the effective critical 
   momentum, thus artificially limiting the avalanche growth rate.
   

Example
-------
Kinetic avalanche runaways can be included in the simulation thusly:

.. code-block:: python

   import DREAM.Settings.Equations.RunawayElectrons as Runaways

   ds = DREAMSettings()
   ...
   ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_KINETIC, pCutAvalanche=0.1)

This will add the Rosenbluth-Putvinski source term to all kinetic equations
in the simulation.

.. note::
   When activating kinetic avalanche, the source integrated 
   over momenta outside the hottail grid, i.e. for all :math:`p>p_\mathrm{max}`, 
   will be added directly to ``n_re``. This way, the total number of knock-ons 
   created is independent of grid parameters.

Synchrotron losses
******************
Synchrotron radiation provides one of the most important energy loss mechanisms
for runaway electrons. DREAM implements the model of
`Decker et al (2016) <https://doi.org/10.1088/0741-3335/58/2/025016>`_, which
uses the radiation time scale

.. math::

   \tau_{\rm rad}^{-1} = \frac{e^4B^2}{6\pi\epsilon_0 m_e^3c^3}.

Synchrotron losses are disabled by default.

List of options
---------------

+------------------------------+-----------------------------------------------------+
| Name                         | Description                                         |
+==============================+=====================================================+
| ``SYNCHROTRON_MODE_INCLUDE`` | Include synchrotron losses in simulation            |
+------------------------------+-----------------------------------------------------+
| ``SYNCHROTRON_MODE_NEGLECT`` | Do **not** include synchrotron losses in simulation |
+------------------------------+-----------------------------------------------------+

Example
-------
Synchrotron losses can be enabled via a call to ``setSynchrotronMode()``:

.. code-block:: python

   import DREAM.Settings.Equations.DistributionFunction as DistFunc

   ds = DREAMSettings()
   ...
   ds.eqsys.f_re.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_INCLUDE)

   # To explicitly disable, call
   #ds.eqsys.f_re.setSynchrotronMode(DistFunc.SYNCHROTRON_MODE_NEGLECT)

Alternatively, a bool may be used to enable/disable synchrotron losses:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.eqsys.f_re.setSynchrotronMode(True)

   # To explicitly disable, call
   #ds.eqsys.f_re.setSynchrotronMode(False)


Momentum threshold
******************

When resolving the thermal population on the grid, the question arises how
to split the electron population into ``cold`` (which enter into the temperature evolution)
and ``hot`` electrons (corresponding to the remainder). This is relevant to the
split of current into ``j_ohm`` and ``j_hot``, density into ``n_cold`` and ``n_hot``,
and also in kinetic equation terms where the contribution from hot electrons are added
as an integral moment of the distribution.

In DREAM, this is handled using the ``pThreshold`` parameter, which can be set using 
three different modes:

+--------------------------------------+--------------------------------------------------------------------------------------------+
| Name                                 | Definition of hot electrons                                                                |
+======================================+============================================================================================+
| ``HOT_REGION_P_MODE_MC``             | :math:`\int_\mathrm{pCutoff}^\mathrm{pMax} \mathrm{d}p`                                    |
+--------------------------------------+--------------------------------------------------------------------------------------------+
| ``HOT_REGION_P_MODE_THERMAL``        | :math:`\int_{\mathrm{pCutoff}\sqrt{2T_\mathrm{cold}/m_e c^2}}^\mathrm{pMax} \mathrm{d}p`   |
+--------------------------------------+--------------------------------------------------------------------------------------------+
| ``HOT_REGION_P_MODE_THERMAL_SMOOTH`` | :math:`\int_0^\mathrm{pMax} \mathrm{d}p ~h(p,\,T_\mathrm{cold})`                           |
+--------------------------------------+--------------------------------------------------------------------------------------------+

.. note::
   The default settings are ``HOT_REGION_P_MODE_THERMAL`` with ``pCutoff = 10``, corresponding to a sharp cutoff
   at ten thermal momenta.

All modes assume an isotropic definition of the hot electrons, where the hot region is independent of angle.
Mode ``MC`` corresponds to a sharp cutoff with ``pThreshold`` the momentum cutoff in units of :math:`m_e c`,
and ``THERMAL`` gives it in units of the thermal momentum. In ``SMOOTH`` mode, the envelope function is defined as:

.. math::
   
   h(p,\,T_\mathrm{cold}) = \tanh\left[\left(p-\mathrm{pCutoff}\sqrt{2T_\mathrm{cold}/m_e c^2}\right)/\Delta p\right]

This smooth envelope function allows the jacobian of the momentum region with respect to the cold electron temperature
to be evaluated, which may improve the rate of convergence in the non-linear solver.

.. warning::
   The step width :math:`\Delta p` in ``SMOOTH`` mode is currently defined as the momentum width of the cell containing the cutoff.
   It is due for a rework, as a failure mode presently occurs when using non-uniform momentum grids with a large jump in step width
   near to the thermal population, in which case a significant fraction of thermal electrons will be counted as `hot`.

Example
-------
The following example shows how to change hot-electron definition to a smooth envelope function located at an elevated cutoff momentum:

.. code-block:: python

   import DREAM.Settings.Equations.HotElectronDistribution as FHot

   ds = DREAMSettings()
   ds.eqsys.f_hot.setHotRegionThreshold(pThreshold=20, pMode=FHot.HOT_REGION_P_MODE_THERMAL_SMOOTH)

This could for example be employed in an attempt to stabilize a hypothetical unstable self-consistent
simulation with a rapid temperature drop when the thermal Maxwellian is resolved on-grid.


Particle source
***************

When running with `COLLISION_FREQUENCY_MODE_FULL`, the entire electron population is followed on the
kinetic grid. A particle source term :math:`S = \mathrm{S\_particle} * S_0(p)` is added to the `hottail` kinetic equation which 
ensures that the correct density evolution is followed; although the kinetic equation is locally density 
conserving, such a source term is needed in order to compensate for changes in the ion charge states by 
ionization and recombination, and also for the radial transport term which is typically constructed such
that it causes a net flow of charge which must be compensated for. 

Currently, the source function shape :math:`S_0(p)` is chosen as the Maxwell-Jüttner distribution evaluated
at the instantaneous cold-electron temperature. It is normalized such that its density moment equals unity, 
so that the unknown quantity `S_particle` represents the net rate at which density is added to (or removed from) 
the distribution.

The particle source is controlled by the following settings:

+------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| Name                         | Model for `S_particle`                                                                                                                      |
+==============================+=============================================================================================================================================+
| ``PARTICLE_SOURCE_ZERO``     | :math:`\mathrm{S\_particle} = 0` (the particle source is neglected)                                                                         |
+------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| ``PARTICLE_SOURCE_IMPLICIT`` | :math:`\int f_\mathrm{hot} \mathrm{d}p = n_\mathrm{cold} + n_\mathrm{hot}` (`S_particle` is such that f_hot has the correct density moment) |
+------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| ``PARTICLE_SOURCE_EXPLICIT`` | :math:`\mathrm{S\_particle} = \Big(\mathrm{d}n_\mathrm{free}/\mathrm{d}t\Big)_\mathrm{ions} + \text{[other sources of density]}`                                    |
+------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------+

Here, `ZERO` deactivates the particle source and allows the density of the hot distribution function
to deviate from `n_cold` + `n_hot`. `IMPLICIT` and `EXPLICIT` yield exactly the same behaviour, but due
to the different formulations yield different condition numbers for the matrix. 

Example
-------

The model used for the particle source can be set by:

.. code-block:: python 

   import DREAM.Settings.Equations.HotElectronDistribution as FHot

   ds = DREAMSettings()

   ds.eqsys.f_hot.setParticleSource(particleSource=FHot.PARTICLE_SOURCE_IMPLICIT):

.. note::
   If the equation system becomes ill-conditioned when running a self-consistent simulation 
   with `COLLFREQ_MODE_FULL`, try changing particle source mode, which may make the equation
   system more well-conditioned.


Tips
****
.. note::

   Are you experiencing unstable solutions? Try to apply one of the
   alternative :ref:`advection term interpolation<advection-interpolation>`
   schemes available in DREAM. These can efficiently suppress numerical
   oscillations and stabilize solutions.
   If the non-linear solver fails to converge, increasing resolution 
   sometimes helps, in particular decreasing the time step.
   If a floating-point exception is thrown, the matrix may have become
   singular to working precision, which is often resolved by using a 
   different linear solver, or increasing grid and/or time resolution.



Class documentation
*******************

.. autoclass:: DREAM.Settings.Equations.DistributionFunction.DistributionFunction
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

.. autoclass:: DREAM.Settings.Equations.HotElectronDistribution.HotElectronDistribution
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

.. autoclass:: DREAM.Settings.Equations.RunawayElectronDistribution.RunawayElectronDistribution
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

