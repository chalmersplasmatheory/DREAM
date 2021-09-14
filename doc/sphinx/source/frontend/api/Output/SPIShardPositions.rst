.. _do-spishardpositions:

SPIShardPositions
=================

The ``SPIShardPositions`` class encapsulates the output of the SPI shard positions :math:`\boldsymbol(x)_\mathrm{p}` obtained from DREAM, and is derived from the :ref:`do-scalarquantity` class. It provides functionality to calculate the radial coordinates from the cartesian coordinates used for the SPI shard positions in DREAM.

.. note::

   This class does currently not have access to the geometry used in the simulation, and hence the calculations made in this class all assume cylindrical geometry.
   
Class documentation
-------------------

.. autoclass:: DREAM.Output.SPIShardPositions.SPIShardPositions
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

