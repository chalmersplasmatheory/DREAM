.. _DREAMSettings:

DREAMSettings
=============

.. toctree::
   :hidden:

   collisions
   eqsys
   momentumgrid
   other
   output
   radialgrid
   solver
   timestep

The ``DREAMSettings`` object represents the settings which are given to a DREAM
simulation as input. It contains a number of convenience and helper routines, as
well as type and option checking, which ensure that many errors are caught even
before DREAM is ever run.


.. contents:: Page overview
   :local:
   :depth: 3


DREAM settings are sorted into a number of categories. Each category contains
settings for a particular module of the kernel:

+---------------------------------------+------------------------------------------------------------------------------+
| Settings category                     | Description                                                                  |
+=======================================+==============================================================================+
| :ref:`ds-collisions`                  | Flags specifying how particle collisions are handled.                        |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`eqsys<ds-eqsys>`                | Settings for all unknown quantities being solved for.                        |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`hottailgrid<ds-momentumgrid>`   | Hot-tail momentum grid settings.                                             |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`other<ds-other>`                | Settings for all "other" quantities (which are not explicitly solved for).   |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`output<ds-output>`              | Settings related to the simulation output.                                   |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`radialgrid<ds-radialgrid>`      | Radial grid settings.                                                        |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`runawaygrid<ds-momentumgrid>`   | Runaway momentum grid settings.                                              |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`solver<ds-solver>`              | Settings for the equation system solver.                                     |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`timestep<ds-timestep>`          | Settings for the module determining the size of time steps to take.          |
+---------------------------------------+------------------------------------------------------------------------------+


Basic usage
-----------
The example below illustrates the basic idea of how to use the ``DREAMSettings``
object for generating a settings file that can be given to the ``dreami``
executable. After running this script, a file named ``dream_settings.h5`` should
appear. By running ``dreami dream_settings.h5``, the simulation can then be
executed.

.. code-block:: python

   from DREAM.DREAMSettings import DREAMSettings
   import DREAM.Settings.CollisionHandler as Collisions

   # Construct a settings object with defaults
   ds = DREAMSettings()

   # Modify settings according to simulation
   ds.collisions.collfreq_mode = Collisions.COLLFREQ_MODE_FULL
   ds.eqsys.E_field.setPrescribedData(0.3)
   ...

   # Save settings to file (which can then be given
   # to the 'dreami' executable)
   ds.save('dream_settings.h5')

Details about which settings are available can be found on the separate pages
linked to in the table at the top of this page.

Chaining simulations
--------------------
Simulations in DREAM can straightforwardly be run in sequence. For example, one
often wants a simulation to start from an initial state which is the solution of
a (possibly non-linear and time-dependent) system of equations. The easiest way
to achieve this may then be to solve the system of equations with DREAM and copy
the solution as input to a subsequent DREAM simulation. An equivalent situation
is when one would like to extend the duration of previous simulation.

To accomodate this need, DREAM simulations can be *chained*. This means that the
settings of the second simulation explicitly refer to the output file of the
first simulation, and load initial values for all quantities from that file.
Specifically, to chain two simulations, one can simply derive the
``DREAMSettings`` object for the second simulation from that of the first:

.. code:: python

   ds1 = DREAMSettings()
   ...
   ds1.output.setFilename('output1.h5')
   ...
   # Chain the two simulations: ds2 takes settings from ds1, and loads
   # initial values from the output of ds1.
   ds2 = DREAMSettings(ds1)

The last line will create a new ``DREAMSettings`` object and transfer all
settings of ``ds1`` into ``ds2``. It will also establish a link between ``ds2``
and the output file ``output1.h5``. This can also be done explicitly, even if
``ds2`` was not created from ``ds1``, using

.. code:: python

   ds2.fromOutput('output1.h5')

Optionally, one can also give ``timeindex`` as the index in the output file of
the time step to take initial values for all quantities from.

Ignoring some initial values
****************************
By default, when two simulations are chained, initial value for all quantities
in the second simulation will be loaded from the output of the first simulation.
This however means that any initial values specified by the user via the
``DREAMSettings`` object will be overridden. To nevertheless initialize one or
more unknown quantity according to the settings (rather than loading their
values from the previous output), one can provide their name in the *ignore
list*:

.. code:: python

   ignorelist = ['N_i', 'W_i']
   ds2.fromOutput('output1.h5', ignore=ignorelist)

Alternatively, can overwrite just the ignore list:

.. code:: python

   ignorelist = ['N_i', 'W_i']
   ds2.setIgnore(ignorelist)

If you want to keep the list of quantities to ignore when chaining simulations,
you need to explicit enforce this when creating the new object:

.. code:: python

   ds2 = DREAMSettings(ds1, keepignore=True)

Copying settings without chaining simulations
*********************************************
Sometimes one may want to copy settings from one ``DREAMSettings`` to another.
This can be achieved in a similar way to when simulations are to be chained:

.. code:: python

   ds2 = DREAMSettings(ds1, chain=False)

This will copy all settings, but will **not** instruct DREAM to take initial
values from the previous output file.

