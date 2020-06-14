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

Initialization
--------------
Initialization of the equation system is handled by the static class
``SimulationGenerator`` (see :ref:`Settings`), which creates the
``EquationSystem`` to use during the simulation, constructs the unknown
quantities and defines their equations.

Since DREAM solves an initial value problem, all unknown quantities must also
be assigned an initial value before beginning the time stepping. Due to the
complicated structure of the DREAM system of equations, where unknown quantities
often depend directly and indirectly on the values of other unknown quantities,
the calculation of unknown quantity initial values is handled by a separate
helper class called ``EqsysInitializer``. The ``EqsysInitializer`` is accessible
from the ``EquationSystem`` as a public property ``EquationSystem::initializer``.

The ``EqsysInitializer`` keeps a list of initialization rules which are run
after the ``EquationSystem`` has been completely constructed, but before the
time stepping begins. Each initialization rule specifies how to initialize the
quantity, as well as which quantities must have been initialized before trying
to set the initial value of the unknown quantity to which the rule applies. The
``EqsysInitializer`` will then calculate an initialization order of the equation
system and execute each initialization rule in order.

DREAM currently allows for two different kinds of initialization:

EVAL_EQUATION
+++++++++++++
This method initializes the unknown quantity by simply evaluating the equation
governing the evolution of the unknown quantity. This requires that the equation
for the unknown is of the form

.. math::

   f(x;y) + g(y) = 0,

where :math:`x` denotes the unknown quantity, :math:`y` is any *other* unknown
quantity (or set of other unknown quantities), :math:`f(x;y)` is an invertible
function of at least :math:`x`, but potentially also other unknown quantities,
and :math:`g(y)` is an arbitrary function of :math:`y` (and importantly NOT of
:math:`x`). An initialization rule of this type will then initialize :math:`x`
by inverting the equation:

.. math::

   x = f^{-1}(-g(y)),

where :math:`f^{-1}(z)` denotes the inverse of `f(x;y)` with respect to
:math:`x`.

EVAL_FUNCTION
+++++++++++++
This method allows the programmer to specify an arbitrary C++ function which
will be used to initialize the unknown quantity. The function takes an
``UnknownQuantityHandler`` as input, as well as a pointer to the value of the
unknown quantity.

STEADY_STATE_SOLVE
++++++++++++++++++
.. warning::

   The initialization method ``STEADY_STATE_SOLVE`` is not yet supported and
   will lead to a fatal error.

   In the future, this method will allow an unknown quantity to be initialized
   by solving the (potentially non-linear) equation governing the evolution of
   the quantity, but with any transient term set to zero.

