ConvergenceScan
===============
The ``ConvergenceScan`` class provides a simple way of quickly generating
parameter scans, particularly in the resolution parameters. A
``ConvergenceScan`` object holds a reference to a baseline :ref:`DREAMSettings`
object, as well as a list of input and output parameters. Each input and output
parameter has a name associate with it, as well as a function which modifies a
given :ref:`DREAMSettings` object (input) or extracts data from a
:ref:`DREAMOutput` object (output).

Once the ``ConvergenceScan`` has been successfully executed, the result can be
conveniently plotted using the :ref:`ConvergenceScanPlot` class.

The most basic use case is illustrated by the following example (although far
more complicated things can be achieved with the class):

.. code-block:: python

   import DREAM
   # ... or
   # from DREAM.ConvergenceScan import ConvergenceScan
   # ... and use as just "ConvergenceScan()"

   # First, create a DREAMSettings object as you normally would...
   ds = DREAM.DREAMSettings()

   # set up ds...
   ...

   # Create convergence scan
   cs = DREAM.ConvergenceScan(settings=ds,
        inparams=['nt', 'hottail.np', 'hottail.nxi'],
        outparams=['other.fluid.runawayRate'])

   # Run convergence scan 
   cs.run()

   # Save results to file
   cs.save('convergence.h5')


.. _convergencescan-by-name:

Set/extract parameter by name
-----------------------------
**TODO**

Scanning until convergence
--------------------------
**TODO**

Object documentation
--------------------

.. autoclass:: DREAM.ConvergenceScan
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__
