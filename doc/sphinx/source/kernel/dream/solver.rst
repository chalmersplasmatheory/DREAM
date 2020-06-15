Solvers
=======
The advance the equation system by one step in time, DREAM can use either a
non-linear solver or solve a linearized form of the equation system using an
Euler-backward time discretization. The non-linear solver is implemented as an
interface to the powerful `SNES <https://www.mcs.anl.gov/petsc/>`_ library
(part of PETSc), while the linearized solver uses the linear equation solvers
of PETSc.

The **solver** in DREAM is invoked once every time step to obtain the solution
for the equation system in a specified time step. The solver will re-build all
equations (and the CollisionQuantityHandler's) one or more times and will return
the solution :math:`\boldsymbol{x}(t_{n+1})` in the new time step. Along with
providing the list of equations to solve, one must also provide the length of
the time step to take, :math:`\Delta t`. The solver does not, in principle,
care about the exact value of :math:`\Delta t`, which makes it possible to
(in theory) adaptively advance the solution in time by separately calculating
the time step to take outside of the solver.

Non-linear and linear forms
---------------------------
Generally, the DREAM equation system can be written on the form

.. math::
   :label: eq_nonlinear

   \boldsymbol{F}(t_n, \boldsymbol{x}) = \boldsymbol{0},

in every time step :math:`t_n`, where :math:`\boldsymbol{x}` represents the
unknown quantities---such as the cold electron density, temperature,
distribution function, ion densities, and more---and
:math:`\boldsymbol{F}(\boldsymbol{x})` is a vector that represents the set of
(non-linear) equations to solve. The **non-linear solver** solves this equation
iteratively using Newton's method:

.. math::

   \boldsymbol{x}_{i+1}(t_n) = \boldsymbol{x}_{i}(t_n) -
   \mathsf{J}^{-1}\left( t_n, \boldsymbol{x}_i \right)
   \boldsymbol{F}\left( t_n, \boldsymbol{x}_i \right).

The non-linear solver comes up with new guesses for the solution
:math:`\boldsymbol{x}_i(t_n)`, which means that DREAM must evaluate both the
function vector :math:`\boldsymbol{F}(t_n,\boldsymbol{x}_i)` and its Jacobian
:math:`\mathsf{J}(t_n,\boldsymbol{x}_i)` in each Newton iteration.

The **linearized solver** starts from the full equation system
:eq:`eq_nonlinear` and linearizes it in time:

.. math::

   \boldsymbol{F}(\boldsymbol{x}_{n+1}) \approx S(\boldsymbol{x}_n) +
   M\left( \boldsymbol{x}_n \right) \boldsymbol{x}_{n+1},

where we use the notation :math:`\boldsymbol{x}_n = \boldsymbol{x}(t_n)` (not
to be confused with the index used for Newton iteration indexing with the
non-linear solver). Here, :math:`M(\boldsymbol{x}_n)` is a matrix and
:math:`S(\boldsymbol{x}_n)` is a vector. The solution of this equation system is

.. math::

   \boldsymbol{x}_{n+1} = -M^{-1}(\boldsymbol{x}_n) S(\boldsymbol{x}_n),

which can be calculated using any of the linear solvers of PETSc.

.. note::

   The linear solver can be significantly more performant than the non-linear
   solver, but it may on the other hand be less accurate if a too large
   :math:`\Delta t\equiv t_{n+1}-t_n` is used in the computation.

