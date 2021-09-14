.. _DREAMOutput:

DREAMOutput
===========

.. toctree::
   :maxdepth: 4
   :hidden:

   Output/DistributionFunction
   Output/EquationSystem
   Output/FluidQuantity
   Output/HotElectronDistributionFunction
   Output/IonHandler
   Output/KineticQuantity
   Output/OtherQuantity
   Output/RunawayElectronDistributionFunction
   Output/ScalarQuantity
   Output/SPIShardRadii
   Output/SPIShardPositions
   Output/UnknownQuantity

The DREAM Python interface contains an elaborate set of classes for handling
DREAM output data.

Base classes
------------

Equation system
---------------

+-------------------------------+----------------------------------------+
| Kinetic quantities            | Description                            |
+===============================+========================================+
| :ref:`f_hot<do-hotelectrons>` | Hot electron distribution function     |
+-------------------------------+----------------------------------------+
| :ref:`f_re<do-reelectrons>`   | Runaway electron distribution function |
+-------------------------------+----------------------------------------+

+--------------------------------------+----------------------------------------+
| Fluid quantities                     | Description                            |
+======================================+========================================+
| :ref:`E_field<do-electricfield>`     | Electric field                         |
+--------------------------------------+----------------------------------------+
| :ref:`j_hot<do-currentdensity>`      | Hot electron current density           |
+--------------------------------------+----------------------------------------+
| :ref:`j_ohm<do-currentdensity>`      | Ohmic current density                  |
+--------------------------------------+----------------------------------------+
| :ref:`j_re<do-currentdensity>`       | Runaway electron current density       |
+--------------------------------------+----------------------------------------+
| :ref:`j_tot<do-currentdensity>`      | Total current density                  |
+--------------------------------------+----------------------------------------+
| :ref:`n_cold<do-fluidquantity>`      | Cold electron density                  |
+--------------------------------------+----------------------------------------+
| :ref:`n_hot<do-fluidquantity>`       | Hot electron density                   |
+--------------------------------------+----------------------------------------+
| :ref:`n_i<do-ionhandler>`            | Ion densities                          |
+--------------------------------------+----------------------------------------+
| :ref:`n_re<do-fluidquantity>`        | Runaway electron density               |
+--------------------------------------+----------------------------------------+
| :ref:`n_tot<do-fluidquantity>`       | Total electron density                 |
+--------------------------------------+----------------------------------------+
| :ref:`psi_p<do-fluidquantity>`       | Poloidal flux                          |
+--------------------------------------+----------------------------------------+
| :ref:`T_cold<do-fluidquantity>`      | Cold electron temperature              |
+--------------------------------------+----------------------------------------+
| :ref:`W_cold<do-fluidquantity>`      | Cold electron energy                   |
+--------------------------------------+----------------------------------------+

+--------------------------------------+----------------------------------------+
| Scalar quantities                    | Description                            |
+======================================+========================================+
| :ref:`I_p<do-scalarquantity>`        | Total plasma current                   |
+--------------------------------------+----------------------------------------+
| :ref:`I_wall<do-scalarquantity>`     | Tokamak wall current                   |
+--------------------------------------+----------------------------------------+
| :ref:`psi_edge<do-scalarquantity>`   | Poloidal flux at plasma edge           |
+--------------------------------------+----------------------------------------+
| :ref:`psi_wall<do-scalarquantity>`   | Poloidal flux at tokamak wall          |
+--------------------------------------+----------------------------------------+
| :ref:`V_loop_w<do-scalarquantity>`   | Loop voltage at tokamak wall           |
+--------------------------------------+----------------------------------------+
| :ref:`x_p<do-spishardpositions>`     | SPI shard positions                    |
+--------------------------------------+----------------------------------------+
| :ref:`Y_p<do-spishardradii>`         | SPI shard radii                        |
+--------------------------------------+----------------------------------------+

Other quantities
----------------



Class documentation
-------------------

.. autoclass:: DREAM.DREAMOutput.DREAMOutput
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __contains__, __getitem__,__init__

