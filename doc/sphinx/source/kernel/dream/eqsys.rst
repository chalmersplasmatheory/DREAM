The Equation System
===================
Much of the DREAM computational kernel is centred around the ``EquationSystem``
class, implemented files under ``src/EquationSystem/``, as well as the header
file ``include/DREAM/EquationSystem.hpp``. The purpose of the equation system
is to keep track of the *unknowns* of the equation system that DREAM solves,
both the equations defined for different unknowns as well as their solution
vectors. As such, the ``EquationSystem`` is one of the most central pieces of
DREAM, and the object from which most operations originate.

The ``EquationSystem`` object has a number of tasks:

- Keep track of the unknowns of the equation system and their equations.
- Rebuild the linear operator matrix/function vector and jacobian for the
  equation system.
- Save solution/ion/grid data at the end of a simulation.

Usage
-----
The following illustrates how an ``EquationSystem`` object is used in DREAM:

0. Construct an ``EquationSystem`` object. The constructor for the
   ``EquationSystem`` takes only the three DREAM grids as input (the fluid
   grid, the hot-tail grid and the runaway grid), as well as the types of the
   momentum grids (:math:`p/\xi` or :math:`p_\parallel/p_\perp`).
1. Populate the equation system.

  - Define unknowns through calls to ``SetUnknown(string name, Grid *grid)``.
  - Define time stepper.

3. Construct and define equations for the various unknowns defined in step 2.
   This is done by successive calls to ``SetEquation()``.
4. Call ``ProcessSystem()`` on the ``EquationSystem`` object -- this primarily
   checks which of the unknowns are *trivial* (and therefore don't need to
   appear in the equation system; see note below) and which are not.
5. Construct the equation system ``Solver()`` object to use for inverting the
   system.


.. note::

   **Trivial unknown**
   is a term used in DREAM to refer to unknowns which do not appear in the
   matrix representation of the equation system. The prime example of such an
   unknown is any prescribed parameter. Since the full evolution of the
   parameter is known ahead of time, no equation needs to solved for the
   unknown, allowing us to keep it out of the equation system matrix.

   The opposite of a *trivial* unknown is a *non-trivial* unknown, and is
   defined as an unknown quantity which must appear in the matrix representation
   of the equation system.

