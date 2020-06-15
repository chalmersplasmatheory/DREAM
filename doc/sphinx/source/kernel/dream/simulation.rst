The Simulation object
=====================
The ``Simulation`` object is the root object of all DREAM simulations. It is the
interface through which one interacts with DREAM in C++. The life-cycle of the
``Simulation`` object looks something like this:

.. code-block:: c++

   // Construct a simulation object from settings
   DREAM::Simulation *sim =
       DREAM::SimulationGenerator::ProcessSettings(settings);

   // Run the full simulation
   sim->Run();

   // Save results to output file
   sim->Save("output.h5");

As this code snippet suggests, the ``Simulation`` object is mostly intended as
a container for the simulation. Configuration should be handled via a call to
the ``ProcessSettings()`` static method of ``SimulationGenerator``, providing a
``Settings`` object with appropriate settings. Output generation is also handled
by the ``Simulation`` object itself with a single call to ``Save()``.
