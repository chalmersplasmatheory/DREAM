
DREAM interface for IMAS
==========================
In this section you will find the details of the DREAM-IMAS interface implementation. Use case leírás

.. toctree::
   :maxdepth: 4
   :hidden:

The IDS input for a basic DREAM simulation is implemented. It requires the following IDS structures

- core_profiles
- equilibrium
- wall

The physical quantities are taken from the IDS fields:

+------------------------------------------------------+---------------------------------------------------------------+
| IDS field                                            | Physical quantity 	                                       |
+======================================================+===============================================================+
| core_profiles/profiles_1d/e_field_parallel           | Electric field (vector)	                               |
+------------------------------------------------------+---------------------------------------------------------------+
| core_profiles/profiles_1d/electrons/temperature      | Electron temperature (vector)                                 |
+------------------------------------------------------+---------------------------------------------------------------+
| core_profiles/profiles_1d/grid/rho_tor               | Minor radius (from core profiles, vector)                     |
+------------------------------------------------------+---------------------------------------------------------------+
| core_profiles/profiles_1d/grid/psi                   | Poloidal flux (from core profiles, vector)                    |
+------------------------------------------------------+---------------------------------------------------------------+
| core_profiles/profiles_1d/ion/density                | Ion densities (vector)                                        |
+------------------------------------------------------+---------------------------------------------------------------+
| core_profiles/vacuum_toroidal_field/r0               | Major radius (scalar)                                         |
+------------------------------------------------------+---------------------------------------------------------------+
| core_profiles/vacuum_toroidal_field/b0               | Magnetic field strength (scalar)                              |
+------------------------------------------------------+---------------------------------------------------------------+
| equilibrium/time_slice/profiles_1d/rho_tor           | Minor radius (from equilibrium, vector)                       |
+------------------------------------------------------+---------------------------------------------------------------+
| equilibrium/time_slice/global_quantities/ip          | Plasma current (scalar)                                       |
+------------------------------------------------------+---------------------------------------------------------------+
| quilibrium/time_slice/profiles_1d/j_tor              | Ohmic current density, (vector)                               |
+------------------------------------------------------+---------------------------------------------------------------+
| wall/description_2d/vessel/unit/annular/resistivity  | Wall resistivity (scalar)                                     |
+------------------------------------------------------+---------------------------------------------------------------+
| equilibrium/time_slice/profiles_2d/psi               | Poloidal flux (from 2D equilibrium, vector)                   |
+------------------------------------------------------+---------------------------------------------------------------+

The wall radius is a user specified quantity which must be specified when calling the IDS read in function. The wall radius and wall resistivity are used to calculate the inverse wall time for the boundary condition of the electric field equation in case a DREAM settings object is initialized.
The radial coordinate used to initialize the DREAM settings object is the mid-plane outer radius. This is acquired by taking the poloidal flux coordinates from the 2D equilibrium IDS on the outer mid-plane. The core profiles radial grid as a function of core profiles poloidal flux is then interpolated onto the outer mid-plane flux to make sure the correct radial grid is used. The plasma profiles located in the core profiles IDS are then interpolated onto the new radial grid in DREAM.

.. code-block:: python

   def readInIDSSlice(shot, run, tokamak, user=os.getlogin(), time=-999, log=False, setUpDream=True, wall_radius=-999):
   
The ``readInIDSSlice()`` function takes a single time slice from the IDS database and loads in the data as an initial condition for a DREAM simulation. The shot and run numbers of the input IDS has to be specified as well as the user from whose IMAS database the shot file will used. The time can also be specified by the user, otherwise the earliest time point for which all the required IDS fields are filled will be loaded. Logging of the IDS loading can be requested. This creates a log file, informing on which of the required fields are loaded successfully.

The ``setUpDREAM()`` option can be used to initialize a DREAM settings object as an output from the function. The object initializes a DREAM run with the initial conditions for the physical quantities set from the IDS and the temperature and electric field calculated self-consistently. A few other settings might be required to be set independently from the IDS read in function. A simple example on how to use IDS's with DREAM can be found in the examples directory in the DREAM repository (``examples/imas``). This is a Python script which creates a ``dream_settings.h5`` file which can be used to run DREAM with the necessary input. The execution of the example file requires the ``imas`` Python module to be installed. The input generation from IDS's is conducted independently from running DREAM with the settings file, and can be run separately.

The option not to create a DREAM settings object is currently not implemented, the loaded data could be saved in Python lists or dictionaries.

.. code-block:: python

   def readInIDS(shot, run, tokamak, user=os.getlogin(), log=False, setUpDream=True, wall_radius=-999):

The ``readInIDS`` fuction works similarly to the time slice function, with the difference of loading the full IDS fields specified in Table 1, instead of a single time slice. A ``DREAMSettings`` object is initialized with prescribed quantities if ``setUpDream()`` is used. 

.. note::

   A simple example on how to use IDS's with DREAM can be found in the examples folder in the DREAM repository (``examples/imas_prescribed``).
