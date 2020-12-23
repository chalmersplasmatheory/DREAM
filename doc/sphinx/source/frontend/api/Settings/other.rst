.. _ds-other:

Other quantities
================
While evolving the main unknown quantities of the equation system, DREAM also
calculates a number of help quantities such as collision frequencies, Coulomb
logarithms, runaway rates etc., which we refer to as "*other quantities*" (since
they are not unknowns). In contrast to most unknowns, these quantities can be
directly evaluated, assuming that the unknowns have been calculated. Since these
quantities are often useful when analyzing simulation output, DREAM provides
options for storing them as well during the calculation.

Due to the large number of quantities available, some of which may occupy
significant amounts of memory, the user must select before the simulation which
quantities to save to the output and which to discard. By default, no other
quantities are stored.

Available quantities
--------------------
.. note::

   The list below may be incomplete. A complete list of available other
   quantities can also be found by browsing the method
   ``OtherQuantityHandler::DefineQuantities()`` in the file
   ``src/OtherQuantityHandler.cpp``.

.. otherquantitylist::
   :source: ../../src/OtherQuantityHandler.cpp


Groups
------

Time grid
---------


