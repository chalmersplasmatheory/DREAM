.. _ds-eqsys-spi:

SPI
===
DREAM includes support for introducing ions to the plasma through Shattered Pellet Injections (SPI). The user can initialise an arbitrary number of shards with given initial sizes, positions, velocities and compositions. This class deals with the settings for the initialization and evolution of the shard sizes, positions and velocities as well as the distribution of the ablated material around the shards. The compositions are set by adding an ion species for every species the pellet is composed of and connecting this ion species to the relevant shards, see the documentation for the :ref:`Ions<ds-eqsys-ions>` class.

.. contents:: Page overview
   :local:
   :depth: 3
   
SPI Model Settings
------------------
There are five types of modes with various options for the SPI modelling: 

* the ablation mode, which determines the model used for the evolution of the shard radii
* the velocity mode, which determines the motion of the shards
* the deposition mode, which determines the spread of the ablated material around the shards
* the heat absorbtion mode, which determines the absorbtion of heat flowing into the neutral cloud surrounding the shards
* the cloud radius mode, which determines the calculation of the radius of the neutral cloud. This radius is used to determine the cross section for the heat absorbtion (if included), and also determines the length scale of the deposition kernel (if the deposition kernel used has a finit width).

This section describes the various modes available in all these five categories.

Ablation Modes
^^^^^^^^^^^^^^
The available ablation modes are

+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| Name                          | Description                                                                                                                                       |
+===============================+===================================================================================================================================================+
| ``ABLATION_MODE_NEGLECT``     | No ablation (default)                                                                                                                             |
+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``ABLATION_MODE_FLUID_NGS``   | Ablation according to the Neutral Gas Shielding (NGS) formula                                                                                     |
|                               | presented by `P. Parks at TSDW 2017 <https://tsdw.pppl.gov/Talks/2017/Lexar/Wednesday%20Session%201/Parks.pdf>`, calculated using the fluid       |
|                               | density and temperature (``T_cold`` and ``n_cold``)                                                                                               |
+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``ABLATION_MODE_KINETIC_NGS`` | Similar to ``ABLATION_MODE_FLUID_NGS``, but expressed in terms of the                                                                             |
|                               | mean energy and isotropic heat flux of the distribution function.                                                                                 |
|                               | NOTE: the factor converting from the total mean energy and isotropic                                                                              |
|                               | heat flux to the corresponding values in the direction towards the                                                                                |
|                               | shards are taken to be the same as for a Maxwellian, so this mode is                                                                              |
|                               | somewhat approximate.                                                                                                                             |
+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

Velocity Modes
^^^^^^^^^^^^^^
The available velocity modes are

+------------------------------+----------------------------------------------+
| Name                         | Description                                  |
+==============================+==============================================+
| ``VELOCITY_MODE NONE``       | The shards do not move (default)             |
+------------------------------+----------------------------------------------+ 
| ``VELOCITY_MODE_PRESCRIBED`` | Prescribed (and constant) shard velocities   |
+------------------------------+----------------------------------------------+


Deposition Modes
^^^^^^^^^^^^^^^^
The available deposition modes are

+------------------------------------------+------------------------------------------------------------------+
| Name                                     | Description                                                      |
+==========================================+==================================================================+
| ``DEPOSITION_MODE_NEGLECT``              | No deposition (default)                                          |  
+------------------------------------------+------------------------------------------------------------------+
| ``DEPOSITION_MODE_LOCAL``                | Delta function source (averaged over the current time step       |
|                                          | and over the grid cells)                                         |
+------------------------------------------+------------------------------------------------------------------+  
| ``DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE`` | Same as ``DEPOSITION_MODE_LOCAL`` but shifted one flux tube      |
|                                          | behind the shards. This can be used to mimic the effect of a     |
|                                          | finite equilibration time in the sense that the ablated material |
|                                          | does not immediately affect the background plasma at its own     | 
|                                          | location, which weakens the feedback effect between the ablation |
|                                          | and the resulting cooling                                        |   
+------------------------------------------+------------------------------------------------------------------+ 
| ``DEPOSITION_MODE_LOCAL_GAUSSIAN``       | Gaussian deposition kernel in the radial coordinate, with        | 
|                                          | 1/e-length scale equal to the radius of the neutral cloud.       |
|                                          | NOTE: This kernel is not averaged over the time step, so it      | 
|                                          | should be used with caution if the shards travel a longer        |
|                                          | distance than the radius of the neutral cloud during a single    |
|                                          | time step!                                                       |
+------------------------------------------+------------------------------------------------------------------+ 

Heat Absorbtion Modes
^^^^^^^^^^^^^^^^^^^^^
The available heat absorbtion modes are

+---------------------------------------------------+-----------------------------------------------------------------+
| Name                                              | Description                                                     |
+===================================================+=================================================================+
| ``HEAT_ABSORBTION_MODE_NEGLECT``                  | No heat absorbtion, only the energy required to ionize          |
|                                                   | the ablated material to the equilibrium distribution of         |
|                                                   | ionisation states is subtracted from the background plasma      |
+---------------------------------------------------+-----------------------------------------------------------------+
| ``HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS``          | All heat flowing through a disc of radius :math:`r_\mathrm{cld}`| 
|                                                   | is subtracted from the background plasma as a delta function    |
|                                                   | source around each shard, where :math:`r_\mathrm{cld}` is the   |
|                                                   | radius of the neutral cloud (no additional sink for the         |
|                                                   | ionization of the ablated material)                             |
+---------------------------------------------------+-----------------------------------------------------------------+
| ``HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN`` | Same as ``HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS`` but with a     |
|                                                   | Gaussian-shaped sink, with 1/e-length scale equal               |
|                                                   | to :math:`r_\mathrm{cld}`. NOTE: This kernel is not averaged    |
|                                                   | over the time step, so it should be used with caution if the    |
|                                                   | shards travel a longer distance than the radius of the          |
|                                                   | neutral cloud during a single time step!                        |
+---------------------------------------------------+-----------------------------------------------------------------+

.. note::

   If the ablated material is immediately deposited locally at the shard positions, they immediately carry the energy absorbed in the cloud back to the plasma, except for the energy needed for the ionization. In such a situation, ``HEAT_ABSORBTION_MODE_NEGLECT`` would be the most physically appropriate heat absorbtion mode.

Cloud Radius Modes
^^^^^^^^^^^^^^^^^^
The available cloud radius modes are

+-------------------------------------------+-------------------------------------------------------------------------+
| Name                                      | Description                                                             |
+===========================================+=========================================================================+
| ``CLOUD_RADIUS_MODE_NEGLECT``             | Cloud radius not used (default)                                         |
+-------------------------------------------+-------------------------------------------------------------------------+
| ``CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT`` | Constant prescribed cloud radius                                        |
+-------------------------------------------+-------------------------------------------------------------------------+
| ``CLOUD_RADIUS_MODE_SELFCONSISTENT``      | Currently gives simply :math:`r_\mathrm{cld}=10r_\mathrm{p}`, which is  |
|                                           | a very rough approximation                                              |
+-------------------------------------------+-------------------------------------------------------------------------+

Adding shards
-------------
The most general way to add shards to the SPI is to directly provide a vector containing the shard sizes, initial positions and velocities to the ``setInitialData()`` method. The initial positions and velocities are given in cartesian coordinates, with the origin at the magnetic axis and the xy-plane coinciding with the poloidal cross section. The vectors specifying the initial positions and velocities should have the format :math:`\boldsymbol{x}_\mathrm{p}=(x_\mathrm{p,1},y_\mathrm{p,1},z_\mathrm{p,1},x_\mathrm{p2},y_\mathrm{p2},z_\mathrm{p2},...)`. The example below shows how to initialise one shard with radius 1 cm on the horisontal mid-plane 2.15 m from the magnetic axis, traveling at a speed of 200 m/s directly towards the magnetic axis.

.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.eqsys.spi.setInitialData(rp=0.01, xp=np.array([2.15,0,0]), vp=np.array([-200,0,0]))
   
There are, however, a number of helper-functions implemented to more easily set up injections with some standard distributions for the shard parameters. These functions are covered in the rest of this section.

Shard positions
^^^^^^^^^^^^^^^
To set ``nShard`` new shard positions to a singel scattering point ``(x0,y0,z0)``, the ``setShardPositionSinglePoint()`` can be used:

.. code-block:: python

   ds.eqsys.spi.setShardPositionSinglePoint(nShard=nShard, shatterPoint=np.array([x0,y0,z0]), add=True)
   
The last argument ``add`` determines wether a new set of shards should be added to the existing ones (``add=True``, default) or if the shard position vector should be cleared (``add=False``).

Shard velocities
^^^^^^^^^^^^^^^^
To select ``nShard`` new shard velocities uniformly within a magnitude range ``v0-DeltaV,v0+DeltaV`` and directions chosen uniformly over an ``nDim`` dimensional cone of opening angle ``alpha/2`` and axis anti-parallell with the x-axis, the setShardVelocitiesUniform() can be used:

.. code-block:: python

   ds.eqsys.setShardVelocitiesUniform(nShard=nShard, abs_vp_mean=v0, abs_vp_diff=DeltaV, nDim=nDim)
   
If ``nDim=1``, all shards simply move along the x-axis. If ``nDim=2``, the direction is chosen along an arc spanning an angle ``alpha`` and if ``nDim=3`` the direction is chosen over an ordinary 3-dimensional cone with opening angle ``alpha/2``. Similarly to ``setShardPositionSinglePoint()``, there is an additional argument ``add`` which is set to ``True`` by default.

In some cases, such as when making a staggered injection, it is practical to create some of the shards and add the corresponding ion species already from the beginning, but not set the shards into motion until a later restart. In this way, one does not have to re-set the ion species when making a restart, but can simply initialise the ion species from the previous output. This can be done by using the additional parameter ``shards`` in ``setShardVelocitiesUniform()``, which is a slice specifying the indices for the shards whose velocities should be updated. If one, for example, wants to change the velocities of the ``nShard2`` last shards added, this could be done as

.. code-block:: python

   ds.eqsys.setShardVelocitiesUniform(abs_vp_mean=v0, abs_vp_diff=DeltaV, nDim=nDim, shards=slice(-nShard2,None))
   
When ``shards`` is not ``None``, the ``nShard`` parameter is automatically set to the number of shards specified by the ``shards``-parameter, and ``add`` is set to ``False``.




