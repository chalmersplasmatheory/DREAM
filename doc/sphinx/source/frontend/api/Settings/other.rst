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

Storing other quantities
------------------------
To store a set of other quantities during a DREAM simulation, you must specify
the names of the quantities to store when setting up the simulation:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.other.include('fluid/Eceff', 'fluid/Ecfree', 'fluid/Ectot', 'fluid/EDreic')

A list of available quantities can be found further down on this page.

.. note::

   The reason that we require the user to specify the name of the other
   quantities desired is to reduce the memory consumption during simulations.
   In particular the kinetic quantities (``hottail/...`` and ``runaway/...``)
   can occupy a significant amount of memory during long simulations.


Available quantities
--------------------

.. otherquantitylist::
   :source: ../../src/OtherQuantityHandler.cpp


Groups
------
Instead of specifying the name of each quantity separately to the
method ``other.include()``, it is also possible to specify the name of *group*
of quantities. The predefined groups are listed below.

.. otherquantitylist::
   :source: ../../src/OtherQuantityHandler.cpp
   :groups: yes

Time grid
---------
The other quantities use almost exactly the same time grid as regular unknowns,
with the only difference that other quantities are not defined for :math:`t=0`.
This is due to the fact that many quantities are only/can only be calculated
while taking/after the first time step. Consequently, you should remove the
first time point in ``grid.t`` whenever plotting or working with other
quantities.

