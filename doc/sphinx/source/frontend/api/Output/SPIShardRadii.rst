.. _do-spishardradii:

SPISardRadii
============

The ``SPIShardRadii`` class encapsulates the output of the SPI shard radii :math:`Y_\mathrm{p}=r_\mathrm{p}^{5/3}` obtained from DREAM, and is derived from the :ref:`do-scalarquantity` class. While the ordinary ``ScalarQuantity``-functions treats :math:`Y_\mathrm{p}` rather than :math:`r_\mathrm{p}`, this class contains wrappers for using the functionality in ScalarQuantity for :math:`r_\mathrm{p}`, and also for the total volume of all shards. For instance, the shard radii for all shards at all time steps can be extracted as follows:

.. code-block:: python

   from DREAM.DREAMOutput import DREAMOutput
   ...
   
   do=DREAMOutput('output.h5')
   rp=do.eqsys.Y_p.calcRadii()
   
A convenient way to illustrate an SPI simulation is to plot the shards in the poloidal cross section at different time frames, perhaps on top of a contour plot of another ``FluidQuantity``, similar to the cover page of `O. Vallhagens MSc thesis <https://ft.nephy.chalmers.se/files/publications/606ddcbc08804.pdf>`_. This can be done with the ``plotPoloidal()`` method (using functionality from :ref:`do-spishardpositions` to calculate the radial coordinates). Moreover, support is included for making animations based on such plots. For example, an animation of an entire simulation with the background showing a contour plot of the total electron density can be created as follows:

.. code-block:: python

   do.eqsys.Y_p.animatePoloidal(backgroundQuantity=do.eqsys.n_tot)
   
.. note::

   This class does currently not have access to the geometry used in the simulation, and hence the ``FluidQuantity`` background data is plotted on cylindrical flux surfaces regardless of the actual geometry.

Class documentation
-------------------

.. autoclass:: DREAM.Output.SPIShardRadii.SPIShardRadii
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__



