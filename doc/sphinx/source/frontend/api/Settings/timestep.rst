.. _ds-timestep:

TimeStepper
===========
The time stepper module is responsible for advancing the system in time. Two
different time steppers are available, namely a fixed length stepper and a
simple adaptive stepper.

.. contents:: Page overview
   :local:
   :depth: 3


Constant step length
--------------------
The constant time stepper advances the system in time at a fixed rate
:math:`\Delta t` set by the user. The user always needs to specify the maximum
simulation time, :math:`t_{\rm max}`, but can choose to either specify the step
length :math:`\Delta t` or the number of steps to take, :math:`N_t`. The
options are set using the following code in the Python interface:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.timestep.setTmax(1e-2)
   ds.timestep.setDt(5e-4)

   # or, alternatively...
   ds.timestep.setTmax(1e-2)
   ds.timestep.setNt(20)


Reducing number of saved steps
******************************
By default, every time step taken is saved to the output file. For long
simulations with fine time resolution this may however require a significant
amount of memory and disk space. To reduce the amount of space needed, the
constant time stepper can instructed to only store a certain number of time
steps in the final output file. This is achieved using the
``setNumberOfSaveSteps()`` method as follows:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.timestep.setNt(6000)
   ds.timestep.setNumberOfTimeSteps(1000)

This will cause the DREAM output file to only contain every 6th time step
taken by the solver. Note that the solutions are unaffected by this
downsampling---the finer time resolution is still used for evolving the
system in time---and the only effect is that the output file contains fewer
of the time steps actually taken.

Adaptive step length
--------------------
A simple adaptive time stepper has been implemented for DREAM. Since DREAM
uses an Euler backward method for time stepping, the adaptive scheme uses basic
two-point time step refinement to estimate the current error in the solution.
The scheme consists of the following steps:

1. Evaluate the solution :math:`\boldsymbol{x}(t+\Delta t)` from :math:`\boldsymbol{x}(t)`
2. Evaluate the solution :math:`\boldsymbol{x}(t+\Delta t/2)` from :math:`\boldsymbol{x}(t)`,
   and then the solution :math:`\boldsymbol{x}(t+\Delta t)` from :math:`\boldsymbol{x}(t+\Delta t/2)`.
3. Estimate the error in each unknown quantity from the difference between the
   two solutions for :math:`t+\Delta t` obtained in steps 1 and 2.
4. If the errors in all monitored unknown quantities are acceptably low, the
   solution is accepted and the optimal time step length for the next step is
   estimated. If the error in one or more of the monitored quantities is too
   large, however, a new optimal time step is estimated based on the error and
   the scheme is repeated for the same time :math:`t` from step 1 above.

.. warning::

   The adaptive time stepper is currently considered too unstable to be used
   in simulations.

.. todo::

   Provide more details about the adaptive time stepper scheme.


Class documentation
-------------------

.. autoclass:: DREAM.Settings.TimeStepper.TimeStepper
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

