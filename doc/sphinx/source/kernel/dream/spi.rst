SPI
===

DREAM supports modelling of material injections in the form of SPI. While there are helper functions available in the python interface to set up standard injection scenarios by providing only a rather small number of parameters, DREAM requires as an input the initial size, composition, position and velocity of every shard. The subsequent evolution of the shards themselves are described by three quantitie; the positions ``x_p`` and velocities ``v_p`` in cartesian coordinates (stored as a vector with the x, y and z-coordinates of a single shard following after each other) and the shard-radii variable ``Y_p``, with :math:`Y_\mathrm{p}=r_\mathrm{p}^{5/3}`, where :math:`r_\mathrm{p}` is the shard radii. The variable ``Y_p`` is used instead of the shard radii to avoid a singularity in the expression for the time derivative of the shard-radii variable when using the Neutral Gas Shielding (NGS) formula for the ablation. These three quantities can then be used to express the other terms related to the SPI.

The currently implemented SPI-related terms include:  

+---------------------------+----------------------------------------------------------------------------------------+
| Name                      | Description                                                                            |
+===========================+========================================================================================+
| ``SPIAblationTerm``       | Decrease in shard sizes due to ablation                                                |
+---------------------------+----------------------------------------------------------------------------------------+
| ``IonSPIDepositionTerm``  | Density increase of (every charge state of) a given ion specie resulting from ablation |
+---------------------------+----------------------------------------------------------------------------------------+
| ``IonSPIIonizLossTerm``   | Ionization energy loss required for the ablation                                       |
+---------------------------+----------------------------------------------------------------------------------------+
| ``SPIHeatAbsorbtionTerm`` | Removal of thermal energy from the background plasma flowing into the pellet cloud     |
+---------------------------+----------------------------------------------------------------------------------------+
| ``SPITransientTerm``      | Transient term adapted to the case of a scalar quantity with multiples                 |
+---------------------------+----------------------------------------------------------------------------------------+

While some calculations are specific for a single term, many calculations are common for several terms. For this reason, the evaluation of the SPI-related terms (and their jacobians) are divided between the ``SPIHandler`` class, performing calculations which concern several terms, and the implementations of each individual term. The individual terms then make use of the calculations performed in the ``SPIHandler`` class, but also add term-specific calculations. 

The SPIHandler class
--------------------
This class constitutes the core of the SPI implementation in DREAM, containing the implementation of (the alternatives for) the ablation rate and deposition kernel. It calculates the time derivative of ``Y_p`` and the resulting number of particles added to every radial grid point, from every shard and in total, as well as the corresponding jacobians. The currently implemented ablation model is the NGS model by `P. Parks presented at TSDW 2017 <https://tsdw.pppl.gov/Talks/2017/Lexar/Wednesday%20Session%201/Parks.pdf>`_, which is valid for pellets consisting of neon and/or deuterium, The user can choose between calculating the ablation based on the Maxwellian electron density and temperature, or by using the mean energy and isotropic heat flux calculated from the distribution function. The deposition kernels available are a delta function kernel, averaged over the time step, and a Gaussian kernel in the radial coordinate (not averaged over the time step!).

The choice of model for the ablation rate and deposition kernel is specified by the input to the constructor. The input to the constructor also includes the pellet composition, which is used to calculate the composition-dependent factor in the ablation rate (as the pellet composition remains constant, this factor only has to be calculated once). 

This class also keeps track of the flux surface coordinates and corresponding indices of each shard, which are needed to evaluate the contribution of every shard to the density on each grid point. The coordinate transformation from the cartesian coordinates itself is however performed by the ``RadialGridGenerator``.

SPIAblationTerm
---------------
This term, together with the ``SPITransientTerm``, comprises the equation for the shard radii variable ``Y_p``. As the ablation rate and corresponing change in ``Y_p`` (and the corresponding jacobian) is already calculated by the ``SPIHandler``, this term only has to collect this information and store it to the appropriate function vector and jacobian matrix block.

IonSPIDepositionTerm
--------------------
This term accounts for the density increase of a specified ion species due to the ablation. This is done by collecting information about the total density at every radial grid point from the ``SPIHandler``, and then multiply this increase by the fraction of the pellet consisting of the ion species corresponding to this term, and also the fraction of the ablated material deposited to the various charge states. The fraction of the pellet consisting of the ion species concerned by this term is given as input to the constructor. As the ionization of the lowest charge states typically occur over a shorter time scale then we want to resolve at the usually high temperature at the start of the injection, the material is deposited into the equilibrium charge state distribution. The corresponding weight factors are thus also calculated by this term.

IonSPIIonizLossTerm
-------------------
This term accounts for the ionization energy required to reach the equilibrium charge state distribution into which the ablated material is deposited.

SPIHeatAbsorbtionTerm
---------------------

SPITransientTerm
----------------




