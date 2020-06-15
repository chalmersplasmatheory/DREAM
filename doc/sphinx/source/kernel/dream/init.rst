.. _Settings:

Settings
========
.. seealso::

   All initialization logic in DREAM can be found under
   `src/Settings/ <https://github.com/chalmersplasmatheory/DREAM/tree/master/src/Settings>`_.

Initialization of the simulation is handled by the static class
``SimulationGenerator`` in ``libdream``. The ``SimulationGenerator``, in turn,
receives it configuration input as a ``Settings`` object. The
``SimulationGenerator`` subsequently constructs a ``Simulation`` object (and all
of its contents), according to the given specification. Calling ``Run()`` on the
``Simulation`` then launches and runs through the entire simulation.

Settings are given to the ``SimulationGenerator`` as a ``Settings`` object,
which can either be constructed directly in the C++ code, or may be constructed
using the ``SFile`` interface with a call to ``DREAM::SettingsSFile::LoadSettings()``.
As such, ``libdream`` does not require settings to be specified in a particular
file format, and settings need not even be provided in a physical file at all
(instead, they could be constructed in an (e.g. C++ or Python) interface
program).

.. note::

   The **SFile** API is part of `SOFTLib <https://github.com/hoppe93/softlib>`_
   and provides a uniform interface for writing (primarily) HDF5 and Matlab
   files.

Settings definitions
--------------------
Before setting configuration options in a ``Settings`` object, the options must
be defined. This is done automatically by ``libdream`` via a call to
``SimulationGenerator::DefineOptions()``. This static method defines all options
which could possibly be specified to DREAM, along with their default values and
brief descriptions of their functionality. By requiring options to be defined
separately from where they are used in the code allows us to construct a list
of all available options, and allows us to search input files for specific
options.

A setting is defined via a call to ``Settings::DefineSetting()``. The method
takes essentially three arguments:

- Name of setting
- Brief setting description
- Default value

A fourth optional argument ``mandatory`` can also be set and forces the setting
to be specified by the user (i.e. no default value will be set). The default
argument above can consist of between 1-3 actual arguments passed to the method,
depending on the argument data type. At the time of writing, the setting may
have one of the following types:

- ``bool`` (flag --- 1 argument: default value)
- ``int_t`` (scalar integer --- 1 argument: default value)
- ``real_t`` (scalar real --- 1 argument: default value)
- ``int_t*`` (integer vector --- 2 arguments: no. of elements, default value)
- ``int_t*`` (multidimensional integer array --- 3 arguments: no. of dimensions, no. of elements, default value)
- ``real_t*`` (real vector --- 2 arguments: no. of elements, default value)
- ``real_t*`` (multidimensional real array --- 3 arguments: no. of dimensions, no. of elements, default value)
- ``std::string`` (string --- 1 argument: default value)

To ensure that the setting is defined with the desired data type, it is
recommended that the default value is explicitly cast to the desired type.

*Several examples of how to define a setting can be found in the files located
under* ``src/Settings/Equations/``. *We provide only a simple example here for
completeness:*


.. code-block:: c++

   /**
    * Define some options in the given 'Settings' object.
    *
    * s: Settings object to define options for.
    */
   void define_settings(Settings *s) {
       const len_t ndims   = 2;
       const len_t dims[2] = {0};

       s->DefineSetting("mymodule/setting1", "Brief description", (int_t)2);
       s->DefineSetting("momodule/setting2", "Brief description", ndims, dims, (real_t*)nullptr);
   }


.. note::

   By convention, we separate settings for different modules by a forward slash
   (``/``). This is however not only an arbitrary choice, as it directly
   translates to a similar grouping in SFile files (HDF5 and Matlab). This means
   that an option with name ``equationsystem/f_hot/n0`` will be stored in a
   Matlab file as a member of the struct ``f_hot``, which in turn is a member of
   the struct ``equationsystem``, and can be accessed as
   ``equationsystem.f_hot.n0`` in Matlab. The SFile automatically makes the
   translation from slash-separated paths to struct objects.

.. warning::

   **Avoid typos!**
   Define module names as macros and use the macros instead of hard-coding
   strings into the ``DefineSetting()`` methods! (see for example
   the ``DefineOptions_Ions()`` in ``src/Settings/Equations/ions.cpp``)

Setting settings
----------------
.. note::

   This section discusses assigning values by hand in C++. Setting values can
   also be loaded from HDF5/Matlab files via a single call to the static method
   ``LoadSettings()`` of ``DREAM::SettingsSFile``.

Values can be assigned to settings via calls to the ``SetSetting()`` methods of
a ``Settings`` object. For every ``DefineSetting()`` method applying to a
specific data type, there is a corresponding ``SetSetting()`` method for
assigning a value to the setting. Each method takes essentially two arguments:

- Name of setting to set
- Setting value

As for the ``DefineSetting()`` methods, the value can consist of up to three
parameters, depending on the data type.

.. code-block:: c++

   void set_settings(Settings *s) {
       const len_t ndims   = 2;
       const len_t dims[2] = {2,2};

       real_t *val = new real_t[dims[0]*dims[1]];
       for (len_t i = 0; i < dims[0]*dims[1]; i++)
           val[i] = i+1;

       s->SetSetting("mymodule/setting1", (int_t)1);
       s->SetSetting("mymodule/setting2", ndims, dims, (real_t*)val);
   }

Reading setting values
----------------------
Settings values can (of course!) also be read from the ``Settings`` object. This
is achieved using the ``GetXXX()`` methods. For each ``DefineSetting()`` and
``SetSetting()`` method, there exists a corresponding getter. The getters take
at least the name of the setting as input. The array getters additionally take
the two parameters ``len_t nExpectedDims`` and ``len_t dims[]`` as input.
Basically, ``dims`` is an array of integers which will contain the size of the
loaded array on return, and ``nExpectedDims`` indicates the size of ``dims``.
If the requested setting is **not** an array with ``nExpectedDims`` dimensions,
a ``SettingsException`` is thrown.

All ``GetXXX()`` methods also take an optional bool argument called
``markused``. If ``true`` (the default), this argument indicates that the
``used`` flag should be set on the setting. This allows DREAM to identify which
settings are actually used in a simulation, and store those settings
specifically along with the simulation output. This can further help the user
identify if their settings are actually recognized by DREAM or not.

.. code-block:: c++

   void get_settings(Settings *s) {
       bool setting1 = s->GetSetting("mymodule/setting1");

       const len_t nExpectedDims = 2;
       len_t dims[2];
       const real_t *setting2 = s->GetSetting("mymodule/setting2", nExpectedDims, dims);

       // Do something with the loaded values
       ...
   }

