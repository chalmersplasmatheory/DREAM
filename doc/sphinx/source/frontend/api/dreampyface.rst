
.. _dreampyface:

dreampyface
===========
The DREAM Python interface is a Python module which allows you to directly call
certain functions in DREAM C++ library. The interface is designed to allow
real-time monitoring of DREAM simulations using Python.

.. contents:: Page overview
   :local:
   :depth: 3

Compiling
---------
Compilation of the ``dreampyface`` module is disabled by default, due to that
it requires special compilation settings which can slightly slow down
simulations. Enabling compilation of the module is however as easy as running
CMake with the flag ``-DDREAM_BUILD_PYFACE=YES`` (as usual, done from the
``build`` directory):

.. code-block:: bash

   cmake .. -DDREAM_BUILD_PYFACE=YES

The next time ``make`` is run, a subdirectory ``build/dreampyface/cxx`` will
be created with the file ``libdreampy.so`` in it. This file is the actual
low-level Python module, and could in principle be imported in your script as
is. We however strongly recommend that you use the higher-level Python API
located directly under ``dreampyface`` in the DREAM root directory, which wraps
the low-level API and provides a more user-friendly interface.

.. note::

   When enabling compilation of ``dreampyface``, you should expect a slight
   performance drop in your simulations, *regardless of whether you run them
   through the Python module or through the usual C++ interface* ``dreami``.
   The reason for this is that certain optimizations (use of jump instructions
   with relative addresses) are not possible for dynamically linked libraries,
   which ``libdreampy.so`` necessarily is. Since the DREAM and FVM libraries
   will be included in ``libdreampy.so``, they must also be compiled with
   so-called position-independent code. This in turn means that since we do
   not want to generate two sets of libraries, also ``dreami`` will be compiled
   with position-independent code.

Using the module
----------------
The first step of using the module is to import it in your Python script.
Assuming that both the ``dreampyface`` directory and
``build/dreampyface/cxx/libdreampy.so`` is in your Python path, importing the
interface is as simple as

.. code-block:: python

   import dreampyface

If the above-mentioned directory and file are not in your Python path, you will
need to add them before importing ``dreampyface``, for example at the top of
your script using

.. code-block:: python

   import sys

   # Root directory of your DREAM installation
   DREAMPATH = '/path/to/DREAM'
   sys.path.append(DREAMPATH)
   sys.path.append(f"{DREAMPATH}/build/dreampyface/cxx")

Running simulations
*******************
The following script illustrates the two most basic ways of running a DREAM
simulation using the ``dreampyface`` API:

.. code-block:: python

   import dreampyface
   from DREAM import DREAMSettings

   ds = DREAMSettings()
   # Setup simulation...
   ...

   # Run simulation
   do = dreampyface.run(ds)
   # ...or
   sim = dreampyface.Simulation(ds)
   do = sim.run()

These two approaches are essentially equivalent---under the hood,
:py:meth:`dreampyface.libdreampy.run` effectively creates a ``Simulation``
object and calls :py:meth:`dreampyface.Simulation.Simulation.run` on it. The
primary difference between the two is that when manually creating the
``Simulation`` object, you will be able to retrieve information about the
simulation *before* running it.

Getting information
*******************
Whenever you have access to a :py:class:`dreampyface.Simulation.Simulation`
object (which can be after manually creating it, or in a
:ref:`callback function<dreampyface-callbacks>`) you will be able to query DREAM
about details of the simulation. The 
:py:class:`dreampyface.Simulation.Simulation` class has a number of methods for
getting information about e.g. simulation length, progress, unknown quantities
being evolved, solution data etc. The following example illustrates one possible
use case, although there are plenty of ways in which the class can be used. For
a complete and up-to-date list of available methods, check the auto-generated
documentation for :py:class:`dreampyface.Simulation.Simulation`.

.. code-block:: python

   import dreampyface

   ...
   sim = dreampyface.Simulation(ds)

   unknowns = sim.unknowns.getInfo()
   for name, info in unknowns.items():
       print('{:12s}  {:8d}  {}'.format(name, info['nelements'], info['description']))

.. _dreampyface-callbacks:

Callbacks
*********
DREAM can be instructed to call Python functions on certain events in the code.
At the moment, two types of events are supported: whenever a time step is
completed, and whenever a non-linear solver iteration is completed.

Registering a callback function (i.e. a function to be called on either of the
aforementioned events) can be done in two ways: either by calling one of the
functions with names starting with ``register_callback_``, or the ``onXXX``
functions of the :py:class:`dreampyface.Simulation.Simulation` class.

.. code-block:: python

   import dreampyface

   ...

   # Procedural interface
   dreampyface.register_callback_iteration_finished(lambda _ : print('Iteration finished'))
   dreampyface.register_callback_timestep_finished(timestepFinished)

   # Object-oriented interface
   sim = dreampyface.Simulation(ds)

   sim.onIteration(lambda _ : print('Iteration finished'))
   sim.onTimestep(timestepFinished)

   # Callback function for when a time step has been completed...
   # ('ptr' is a pointer to the C++ Simulation object)
   def timestepFinished(ptr):
       s = dreampyface.Simulation(ptr)
       print('Time: {} / {}'.format(s.getCurrentTime(), s.getMaxTime()))

Implementation details
----------------------
The Python bindings for libdream are written using both C++ and Python. The
C++ code uses the CPython C API to provide a number of functions to the Python
code for interacting with a C++ Simulation object (which is the "mother" object
in any DREAM simulation and provides access to every part of the simulation).

The interface contains C++ code for translating a Python dictionary into a
DREAM ``Settings`` object, as well as an ``SFile`` class for reading/writing a
Python dictionary as if it were a data file.

API Reference
-------------

libdreampy
**********

.. automodule:: dreampyface.libdreampy
   :members:

Simulation
**********

.. autoclass:: dreampyface.Simulation.Simulation
   :members:
   :undoc-members:
   :special-members: __init__

UnknownQuantityHandler
**********************

.. autoclass:: dreampyface.UnknownQuantityHandler.UnknownQuantityHandler
   :members:
   :undoc-members:
   :special-members: __init__

