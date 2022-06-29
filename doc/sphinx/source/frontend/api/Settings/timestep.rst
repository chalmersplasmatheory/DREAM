.. _ds-timestep:

TimeStepper
===========
The time stepper module is responsible for advancing the system in time. Three
different time steppers are available, namely a fixed length stepper, a simple
error estimating adaptive stepper, and an adaptive stepper which estimates and
rescales the time step based on the current ionization time.

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
   ds.timestep.setNumberOfSaveSteps(1000)

This will cause the DREAM output file to only contain every 6th time step
taken by the solver. Note that the solutions are unaffected by this
downsampling---the finer time resolution is still used for evolving the
system in time---and the only effect is that the output file contains fewer
of the time steps actually taken.

Error-based adaptive step length
--------------------------------
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

Ionization-based adaptive step length
-------------------------------------
Instabilities in DREAM primarily arise when the non-linear system of equations
is evolved with too large time steps, so that the initial guess of Newton's
method (which is that the solution is the same as in the previous time step) is
too far from the true solution. These instabilities are therefore caused by
physical processes which cause significant variation on time scales shorter than
the resolved time scale.

One of the shortest time scales usually encountered in DREAM is that of the
ionization of neutral and partially ionized atoms. As greater fractions of the
plasma becomes ionized, however, the ionization time scale increases by several
orders of magnitude. Therefore, if a simulation is started with a time step
that respects the early ionization time scale, it is likely that after a while,
the time step is much smaller than necessary, causing the simulation to take
unnecessarily long time.

With the ionization-based adaptive time stepper, we continuously update the time
step length :math:`\Delta t` as the ionization time scale evolves, by estimating
the current ionization time scale from the variation in the cold electron
density:

.. math::

   \frac{1}{\tau_{\rm ioniz}} \approx \frac{1}{n_{\rm cold}^{(k-1)}}\frac{n_{\rm cold}^{(k)}-n_{\rm cold}^{(k-1)}}{\Delta t}

The user then provides the initial time step :math:`\Delta t_0` at :math:`t=0`,
after which the time stepper automatically rescales the time step throughout the
simulation. A maximum time step can also be selected, since eventually other
physical time scales will be shorter than the ionization time scale, making a
constant time step suitable.

The time stepper can be instructed to automatically determine the first time
step. This is achieved by taking a first very short time step, after which the
ionization time scale can be estimated. Since the required time step length will
depend also on the numerics (e.g. a kinetic simulation could require higher time
resolution than a fluid simulation, since the ionization could induces "greater"
variations in the distribution function), a safety factor can be used to adjust
how the initial time step is determined.

Options
*******
To use the ionization-based adaptive time stepper, one should call the
``setIonization()`` method on the time stepper object. If the user manually
specifies the initial time step length, the method takes the following
arguments:

+------------+----------------------------+
| Option     | Description                |
+============+============================+
| ``dt0``    | Initial time step.         |
+------------+----------------------------+
| ``dtmax``  | Maximum time step allowed. |
+------------+----------------------------+
| ``tmax``   | Final time of simulation.  |
+------------+----------------------------+

If the time stepper should automatically determine the initial time step, the
method takes the following arguments (should be given as keyword arguments):

+-------------------+------------------------------+-------------------------------------------------------------------------------------+
| Option            | Default value                | Description                                                                         |
+===================+==============================+=====================================================================================+
| ``automaticstep`` | :math:`10^{-12}\,\mathrm{s}` | Time step used first to estimate the initial ionization time scale.                 |
+-------------------+------------------------------+-------------------------------------------------------------------------------------+
| ``dtmax``         | None                         | Maximum time step allowed.                                                          |
+-------------------+------------------------------+-------------------------------------------------------------------------------------+
| ``safetyfactor``  | 50                           | Factor giving the relation between the ionization time scale and initial time step. |
+-------------------+------------------------------+-------------------------------------------------------------------------------------+
| ``tmax``          | Final time of simulation.    | Final time of simulation.                                                           |
+-------------------+------------------------------+-------------------------------------------------------------------------------------+

Examples
********
To manually specify the initial time step, the following call can be made:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.timestep.setIonization(dt0=1e-10, dtmax=1e-5, tmax=0.003)

Alternatively, ``tmax`` can be set separately:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.timestep.setIonization(dt0=1e-10, dtmax=1e-5)
   ds.timestep.setTmax(0.003)

If the time stepper should automatically determine the initial time step, the
following call can be made:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.timestep.setIonization(dtmax=1e-5, tmax=0.003)

Alternatively, if you want to explicitly set the time step used to determine the
ionization time scale, or if you want to change the safety factor used to relate
the ionization time scale to the initial time step, you can make the call:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.timestep.setIonization(automaticstep=1e-12, safetyfactor=50, dtmax=1e-5, tmax=0.003)

Class documentation
-------------------

.. autoclass:: DREAM.Settings.TimeStepper.TimeStepper
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

