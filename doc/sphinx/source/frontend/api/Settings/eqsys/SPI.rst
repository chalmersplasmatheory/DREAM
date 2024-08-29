.. _ds-eqsys-spi:

SPI
===
DREAM includes support for introducing ions to the plasma through Shattered Pellet Injections (SPI). The user can initialise an arbitrary number of shards with given initial sizes, positions, velocities and compositions. This class deals with the settings for the initialization and evolution of the shard sizes, positions and velocities as well as the distribution of the ablated material around the shards. The compositions are set by adding an ion species for every species the pellet is composed of and connecting this ion species to the relevant shards, see the documentation for the :ref:`Ions<ds-eqsys-ions>` class.

.. contents:: Page overview
   :local:
   :depth: 3
   
SPI Model Settings
------------------
There are six types of modes with various options for the SPI modelling: 

* the :ref:`ablation mode<ds-eqsys-spi-ablation>`, which determines the model used for the evolution of the shard radii
* the :ref:`velocity mode<ds-eqsys-spi-velocity>`, which determines the motion of the shards
* the :ref:`deposition mode<ds-eqsys-spi-deposition>`, which determines the spread of the ablated material around the shards
* the :ref:`shift mode<ds-eqsys-spi-shift>`, which determines the mode used for the shift of the ablated material due to the plasmoid drift effect
* the :ref:`heat absorbtion mode<ds-eqsys-spi_heatabsorbtion>`, which determines the absorbtion of heat flowing into the neutral cloud surrounding the shards
* the :ref:`cloud radius mode<ds-eqsis-spi-cloudradius>`, which determines the calculation of the radius of the neutral cloud. This radius is used to determine the cross section for the heat absorbtion (if included), and also determines the length scale of the deposition kernel (if the deposition kernel used has a finit width).

This section describes the various modes available in all these five categories.

.. _ds-eqsys-spi-ablation:

Ablation Modes
^^^^^^^^^^^^^^
The available ablation modes are

+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| Name                          | Description                                                                                                                                       |
+===============================+===================================================================================================================================================+
| ``ABLATION_MODE_NEGLECT``     | No ablation (default)                                                                                                                             |
+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``ABLATION_MODE_FLUID_NGS``   | Ablation according to the Neutral Gas Shielding (NGS) formula                                                                                     |
|                               | presented by `P. Parks at TSDW 2017 <https://tsdw.pppl.gov/Talks/2017/Lexar/Wednesday%20Session%201/Parks.pdf>`_, calculated using the fluid      |
|                               | density and temperature (``T_cold`` and ``n_cold``)                                                                                               |
+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``ABLATION_MODE_KINETIC_NGS`` | Similar to ``ABLATION_MODE_FLUID_NGS``, but expressed in terms of the                                                                             |
|                               | mean energy and isotropic heat flux of the distribution function.                                                                                 |
|                               | NOTE: the factor converting from the total mean energy and isotropic                                                                              |
|                               | heat flux to the corresponding values in the direction towards the                                                                                |
|                               | shards are taken to be the same as for a Maxwellian, so this mode is                                                                              |
|                               | somewhat approximate.                                                                                                                             |
+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

.. _ds-eqsys-spi-velocity:

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

.. _ds-eqsys-spi-deposition:

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

.. _ds-eqsys-spi-shift:

Shift Modes
^^^^^^^^^^^^^^^^
The available shift modes are

+------------------------------------------+---------------------------------------------------------------------------+
| Name                                     | Description                                                               |
+==========================================+===========================================================================+
| ``SHIFT_MODE_NEGLECT``                   | No shift (default)                                                        |  
+------------------------------------------+---------------------------------------------------------------------------+
| ``SHIFT_MODE_ANALYTICAL``                | Calculate the shift using the analytical model derived by                 |
|                                          | `Vallhagen et al, JPP 2023 <https://doi.org/10.1017/S0022377823000466>`_  |
+------------------------------------------+---------------------------------------------------------------------------+  

.. _ds-eqsys-spi-heatabsorbtion:

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

.. _ds-eqsys-spi-cloudradius:

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

   ds.eqsys.setShardVelocitiesUniform(nShard=nShard, abs_vp_mean=v0, abs_vp_diff=DeltaV, alpha_max=alpha, nDim=nDim)
   
If ``nDim=1``, all shards simply move along the x-axis. If ``nDim=2``, the direction is chosen along an arc spanning an angle ``alpha`` and if ``nDim=3`` the direction is chosen over an ordinary 3-dimensional cone with opening angle ``alpha/2``. Similarly to ``setShardPositionSinglePoint()``, there is an additional argument ``add`` which is set to ``True`` by default.

In some cases, such as when making a staggered injection, it is practical to create some of the shards and add the corresponding ion species already from the beginning, but not set the shards into motion until a later restart. In this way, one does not have to re-set the ion species when making a restart, but can simply initialise the ion species from the previous output. This can be done by using the additional parameter ``shards`` in ``setShardVelocitiesUniform()``, which is a slice specifying the indices for the shards whose velocities should be updated. If one, for example, wants to change the velocities of the ``nShard2`` last shards added, this could be done as

.. code-block:: python

   ds.eqsys.setShardVelocitiesUniform(abs_vp_mean=v0, abs_vp_diff=DeltaV, alpha_max=alpha, nDim=nDim, shards=slice(-nShard2,None))
   
When ``shards`` is not ``None``, the ``nShard`` parameter is automatically set to the number of shards specified by the ``shards``-parameter, and ``add`` is set to ``False``.

Shard radii
^^^^^^^^^^^
There are a number of helper methods implemented to select shard sizes from the Bessel-like statistical distribution found in `P. Parks GA report <10.2172/1344852>`_:

.. math::
   P(r_\mathrm{p}) = k_\mathrm{p}^2 r_\mathrm{p} K_0(k_\mathrm{p} r_\mathrm{p}),

where :math:`K_0` is the zeroth modified Bessel function of the second kind. The inverse characteristic shard size :math:`k_\mathrm{p}` is related to the particle content, composition and degree of shattering of the pellet according to

.. math::
   k_\mathrm{p}=\left(\frac{6\pi^2n_\mathrm{p}N_\mathrm{s}}{N_\mathrm{inj}}\right)^{1/3},

where :math:`n_\mathrm{p}` is the solid particle density of the pellet, :math:`N_\mathrm{s}` is the number of shards into which the pellet is shattered and :math:`N_\mathrm{inj}` is the total number of particles contained in the pellet.

If the desired :math:`k_\mathrm{p}` is already known, one can sample ``nShard`` shards from the above distribution as

.. code-block:: python

   rp=ds.eqsys.spi.sampleRpDistrParksStatistical(N=nShard, kp=kp)

If one instead wants to specify :math:`N_\mathrm{s}`, :math:`N_\mathrm{inj}` and the pellet composition, one can instead use the ``setRpParksStatistical()``-method. The following lines of code illustrate how to add a pellet containing :math:`10^{24}` particles shattered into 1000 shards, consisting of 5% neon and 95% deuterium, with radii selected from the above distribution. The Ion species connected to this pellet are also added by the ``setRpParksStatistical()``-method by passing a pointer to the ``ds.eqsys.n_i``-object.
 
.. code-block:: python

   Ninj=1e24 # Total number of injected particles
   nShard=1000 # Number of shards into which the pellet is shattered
   Zs=[1,10] # List of charge numbers of the species the pellet is composed of
   isotopes=[2,0] # List of isotopes, 0 meaning naturally occuring mix
   molarFractions=[0.95,0.05] # List of molar fractions specifying the pellet composition
   ionNames=['D_inj','Ne_inj'] # List of names of the ion species connected to this pellet

   ds.eqsys.spi.setRpParksStatistical(nShard=nShard, Ninj=Ninj, Zs=Zs, isotopes=isotopes, molarFractions=molarFractions, ionNames=ionNames, n_i=ds.eqsys.n_i)

As earlier, an extra argument ``add=False`` resets the shard radii instead of adding new ones.


Example: staggered injection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The code block below illustrates how to set up staggered Deuterium-Neon injections similar to those investigated in `O. Vallhagens MSc thesis <https://ft.nephy.chalmers.se/files/publications/606ddcbc08804.pdf>`_, using the wrapper function setParamsVallhagenMSc() to set the initial positions, velocities and radii of the shards on the same line.

.. code-block:: python

   radius=[0,2] # Span of the radial grid
   radius_wall=2.15 # Location of the wall
   ...
   
   # Settings for the first deuterium SPI
   nShardD=1742 # Number of shards
   NinjD=2e24 # Number of atoms
   alpha_maxD=0.17 # Divergence angle
   abs_vp_meanD=800 # Mean shard speed
   abs_vp_diffD=0.2*abs_vp_meanD # Width of the uniform shard speed distribution

   ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardD, Ninj=NinjD, Zs=[1], isotopes=[2], molarFractions=[1], ionNames=['D_inj'], n_i=ds.eqsys.n_i, shatterPoint=np.array([0,0,radius_wall]), abs_vp_mean=abs_vp_meanD, abs_vp_diff=abs_vp_diffD, alpha_max=alpha_maxD)

   # Settings for the second neon SPI
   nShardNe=50 # Number of shards
   NinjNe=1e23 # Number of atoms
   alpha_maxNe=0.17 # Divergence angle
   abs_vp_meanNe=200 # Mean shard speed
   abs_vp_diffNe=0.2*abs_vp_meanNe # Width of the uniform shard speed distribution

   # Initialise neon shards with zero velocity, and set the velocity before the restart when the neon injection should actually start
   ds.eqsys.spi.setParamsVallhagenMSc(nShard=nShardNe, Ninj=NinjNe, Zs=[10], isotopes=[0], molarFractions=[1], ionNames=['Ne_inj'], n_i=ds.eqsys.n_i, shatterPoint=np.array([0,0,radius_wall]), abs_vp_mean=0, abs_vp_diff=0, alpha_max=alpha_maxNe)

   ...
   runiface(ds, 'output_D_inj.h5')
   ...
   
   ds2=DREAMSettings(ds)
   ds2.fromOutput('output_D_inj.h5', ignore=['v_p','x_p'])
   ...
   
   # Set the shards of the second injection into motion and advance them until
   # the fastest shards reach the plasma edge
   do=DREAMOutput('output_D_inj.h5')
   ds2.eqsys.spi.vp=do.eqsys.v_p.data[-1,:].flatten()
   ds2.eqsys.spi.xp=do.eqsys.x_p.data[-1,:].flatten()
                
   ds2.eqsys.spi.setShardVelocitiesUniform(abs_vp_mean=abs_vp_meanNe,abs_vp_diff=abs_vp_diffNe,alpha_max=alpha_maxNe,shards=slice(-nShardNe,None))

   t_edge=(radius_wall-radius[1])/np.max(-ds2.eqsys.spi.vp[-3*nShardNe::3])
   ds2.eqsys.spi.xp[-3*nShardNe:]=ds2.eqsys.spi.xp[-3*nShardNe:]+ds2.eqsys.spi.vp[-3*nShardNe:]*t_edge

   runiface(ds2, 'output_Ne_inj.h5')
   
Plasmoid drift
--------------
The effect of plasmoid drifts can be accounted for using the analytical model given by equation (A4) in `Vallhagen et al, JPP 2023 <https://doi.org/10.1017/S0022377823000466>`_. To use this model, the following model parameters have to be speified: representative temperature :math:`T` during the drift duration; the initial temperature :math:`T_0` close to the pellet (before the cloud has drifted away from the pellet); the characteristic half-thickness of the drifting cloud :math:`\Delta y` (should be similar to the radius of the neutral cloud around the pellet); the major radius :math:`R_\mathrm{m}` at the magnetic axis (only used if the major radius is otherwise infinite in the simulation); the average charge states :math:`Z_\mathrm{avg}^{(i)}` inside the drifting cloud for all ion species include in the pellets. Note that the average charge states can not be calculated using the ADAS rates because the conditions in the drifting cloud, especially the density and optical thickness, are very different from the validity range and assumptions in ADAS, and we therefore take user-given estimates for them. These estimates should be given as a look-up table in the form of the lists ``ZavgArray``, ``Zs`` and ``isotopes``, where ``ZavgArray`` provides the average charge states to be used for ion species with atomic numbers ``Zs`` and isotopes ``isotopes``.


Below is an example showing how to include the analytical drift model for a mixed D+Ne SPI with typical settings. The values of the model parameters are motivated by e.g. numerical studies by `Matsuyama <https://doi.org/10.1063/5.0084586>`_ and experimental studies by `Müller et al <https://doi.org/10.1088/0029-5515/42/3/311>`_. In particular, these studies indicate that the representative temperature :math:`T` is significantly different for purely hydrogenic pellets and impurity doped pellets, as the doped pellets will have much stronger radiation losses. Therefore, different temperatures can be set
for different shards as long as they are in the correct order. The following example includes ``nShardD``
deuterium shards (atomic number 1, isotope 2)  with :math:`T=30\,\rm eV` and ``nShardNe`` number of neon doped shards (atomic number 10, isotope 0 i.e. naturally occuring mix) with :math:`T=5\,\rm eV`.

.. code-block:: python

    ds.eqsys.spi.setShiftParamsAnalytical(shift = SHIFT_MODE_ANALYTICAL, T = np.concatenate([30*np.ones(nShardD), 5*np.ones(nShardNe)]), T0 = 2, delta_y = 0.0125, Rm = R0 , ZavgArray =[1 , 2], Zs, [1,10], isotopes = [2,0]):

Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.SPI.SPI
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

