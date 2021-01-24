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
   solver

The ``DREAMSettings`` object represents the settings which are given to a DREAM
simulation as input. It contains a number of convenience and helper routines, as
well as type and option checking, which ensure that many errors are caught even
before DREAM is ever run.

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
| ``radialgrid``                        | Radial grid settings.                                                        |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`runawaygrid<ds-momentumgrid>`   | Runaway momentum grid settings.                                              |
+---------------------------------------+------------------------------------------------------------------------------+
| :ref:`solver<ds-solver>`              | Settings for the equation system solver.                                     |
+---------------------------------------+------------------------------------------------------------------------------+
| ``timestep``                          | Settings for the module determining the size of time steps to take.          |
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
