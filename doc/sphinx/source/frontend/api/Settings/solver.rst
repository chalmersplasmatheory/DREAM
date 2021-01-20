.. _ds-solver:

Solver
======
The *solver* is the part of the code that is responsible for solving a given
system of non-linear equations, and thereby computing the system state in the
next time step. DREAM provides two different solvers, namely a linearly implicit
solver as well as a regular Newton solver. In the discussion that follows, we
let :math:`\boldsymbol{F}(t,\boldsymbol{x(t)}) = 0` represent the non-linear
equation system that must be solved to advance the system in time.

.. contents:: Page overview
   :local:
   :depth: 3

Which solver to use is determined by calling ``setType()`` on the solver object:

.. code-block::

   import DREAM.Settings.Solver as Solver

   ds = DREAMSettings()
   ...
   ds.solver.setType(Solver.NONLINEAR)

The available solver types are

+---------------------+------------------------------+
| Name                | Description                  |
+=====================+==============================+
| ``LINEAR_IMPLICIT`` | The linearly implicit solver |
+---------------------+------------------------------+
| ``NONLINEAR``       | The Newton solver            |
+---------------------+------------------------------+

Linearly implicit solver
------------------------
The linearly implicit solver linearizes the equation system in time around the
next time, :math:`t_{l+1}`, as

.. math::

   \boldsymbol{F}\left( t_{l+1}, \boldsymbol{x}^{(l+1)} \right) \approx
   \mathsf{M}\left( t_l, \boldsymbol{x}^{(l)} \right) \boldsymbol{x}^{(l+1)} +
   \boldsymbol{S}\left( t_l, \boldsymbol{x}^{(l)} \right),

where :math:`\boldsymbol{x}^{(l)}\equiv\boldsymbol{x}(t_l)`, :math:`\mathsf{M}`
is a matrix and :math:`\boldsymbol{S}` is a vector. Setting both sides to zero
allows us to obtain the solution in the next time step by inverting the matrix
:math:`\mathsf{M}`:

.. math::

   \boldsymbol{x}^{(l+1)} = \mathsf{M}^{-1}\left( \boldsymbol{x}^{(l)} \right)
   \boldsymbol{S}\left( \boldsymbol{x}^{(l)} \right)


Non-linear solver
-----------------
The non-linear solver linearizes the equation system around a point
:math:`\tilde{\boldsymbol{x}}` in the vicinity of the true solution (rather than
around a particular time :math:`t`, as the linearly implicit solver), yielding

.. math::

   \boldsymbol{F}\left(\boldsymbol{x}\right) \approx
   \boldsymbol{F}\left(\tilde{\boldsymbol{x}}\right) +
   \mathsf{J}(\tilde{\boldsymbol{x}})\left( \boldsymbol{x} - \tilde{\boldsymbol{x}} \right)

Setting :math:`\boldsymbol{F}(\boldsymbol{x})` and solving for
:math:`\boldsymbol{x}` then yields the iterative *Newton's method*:

.. math::

   \boldsymbol{x}^{(l+1)}_{i+1} = x^{(l+1)}_i - \mathsf{J}^{-1}\left(
   \boldsymbol{x}^{(l+1)}_i \right)
   \boldsymbol{F}\left( \boldsymbol{x}^{(l+1)}_i \right),

where :math:`\boldsymbol{x}^{(l+1)}_i` denotes the :math:`i` th approximation
to the solution :math:`\boldsymbol{x}(t_{l+1})`.

Which solver should I use?
--------------------------
The difference between the two solvers is primarily that the linearly implicit
solver is fast, but potentially inaccurate and prone to crashing if the time
evolution is rapid or highly non-linear, while the non-linear solver is
generally slower but very robust and guarantees a certain accuracy in the
solution. In principle it should be possible to use either solver to solve any
given problem, although certain non-linear systems have proven difficult for
the linearly implicit solver.

.. tip::

   When in doubt, we recommend using the Newton (non-linear) solver. While it is
   generally slower than the linearly implicit solver, we often find that its
   robustness makes up for this.

In practice, one should therefore consider a number of questions when deciding
which solver to use:

- **Is the equation system non-linear?** If it is, the Newton solver will often
  be required to avoid numerical instabilities crashing the linearly implicit
  solver.
- **Does the system involve rapidly time-varying parameters?** If so, the
  linearly implicit solver may require exceedingly short time steps which could
  partly be avoided by the more accurate Newton solver.


Linear solvers
--------------
Both the linearly implicit and non-linear solvers involve solutions of sets of
linear equations. Hence, both solvers utilize third-party linear equation
solvers available through the PETSc library. Both solvers currently support the
use of four different LU factorization algorithms, namely

+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Name                      | Description                                                                                                                                                                                                                                                                       |
+===========================+===================================================================================================================================================================================================================================================================================+
| ``LINEAR_SOLVER_LU``      | The LU direct solver available in `PETSc <https://www.mcs.anl.gov/petsc/>`_.                                                                                                                                                                                                      |
+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``LINEAR_SOLVER_MKL``     | The parallel direct solver `PARDISO <https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-interface.html>`_, part of Intel's Math Kernel Library (MKL).  |
+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``LINEAR_SOLVER_MUMPS``   | The `MUMPS <http://mumps.enseeiht.fr/>`_ parallel sparse direct solver.                                                                                                                                                                                                           |
+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``LINEAR_SOLVER_SUPERLU`` | The `SuperLU <https://portal.nersc.gov/project/sparse/superlu/>`_ direct LU solver.                                                                                                                                                                                               |
+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
