.. _ds-eqsys:

EquationSystem
==============
The ``EquationSystem`` class is keeps track of settings for all unknown
quantities which DREAM solves for. This class does very little work on its own
and rather acts as a container for the objects representing specific unknown
quantities.

Each unknown quantity which has specific settings is represented by its own
Python class. The object is initialized with defaults when the
``EquationSystem`` class is created (which happens as soon as the
parent ``DREAMSettings`` object is created) and can be accessed as a property
of the equation system. The following unknown quantities are available in any
``EquationSystem`` object:

.. toctree::
   :hidden:

   eqsys/efield

+------------------------------------+------------------------------------+
| Quantity                           | Description                        |
+====================================+====================================+
| :ref:`E_field<ds-eqsys-E_field>`   | Electric field                     |
+------------------------------------+------------------------------------+
| :ref:`f_hot<ds-eqsys-f_hot>`       | Hot electron distribution function |
+------------------------------------+------------------------------------+
| :ref:`n_cold<ds-eqsys-n_cold>`     | Cold electron density              |
+------------------------------------+------------------------------------+
| :ref:`n_i<ds-eqsys-ions>`          | Ion densities and charge states    |
+------------------------------------+------------------------------------+
| :ref:`T_cold<ds-eqsys-T_cold>`     | Cold electron temperature          |
+------------------------------------+------------------------------------+

Examples
--------
Unknown quantities can be accessed in the following way:

.. code-block:: python

   ds = DREAMSettings()

   ds.eqsys.E_field.setPrescribedData(0.3)  # Uniform electric field profile (V/m)
   ds.eqsys.T_cold.setPrescribedData(1100)  # Uniform temperature profile (eV)

   ...

