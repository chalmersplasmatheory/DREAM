.. _ds-static:

Static quantities
=================
Similarly to :ref:`other quantities <ds-other>`, static quantities are
helper quantities which are calculated during the simulation, but are not stored
as regular unknown quantities. In contrast to other quantities, static
quantities do not evolve in time, and so they are fixed throughout the
simulation.

Like other quantities, the user must explicitly indicate which static quantities
to store from a simulation. By default, no static quantities are stored.

Storing static quantities
-------------------------
To store a set of static quantities during a DREAM simulation, you must specify
the names of the quantities to store when setting up the simulation:

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.static.include('boltzmann/deltaAngleKernel_hot_hot')

A list of available quantities can be found further down on this page.

Available quantities
--------------------

.. otherquantitylist::
   :source: ../../src/StaticQuantityHandler.cpp


