General overview
================
The DREAM library, ``libdream``, contains most, if not all, of the physics
available in DREAM.

Key concepts
============
This section describes a number of key concepts within the DREAM library.
The following are concepts which all developers of DREAM should be familiar
with:

.. toctree::
   :maxdepth: 2

   eqsys.rst
   init.rst
   ions.rst
   solver.rst
   spi.rst

Helper classes
==============
This section describes various helper classes which may be useful to
programmers, but which most DREAM programmers will not interact with.

.. toctree::
   :maxdepth: 2

   adas.rst
   constants.rst
   simulation.rst

Equation terms
==============
While the fundamental equation logic is implemented in the FVM library, few of
the equation terms implemented there are directly used in the DREAM library.
Instead, new equation terms are derived from the classes provided in the FVM
library, with specific physics added to the implementations in ``libdream``.
Below is a list of equation terms implemented in DREAM (it aims at being
complete, but is updated manually and may as such be only partial).

.. toctree::
   :maxdepth: 2

   equations/fluid
   equations/kinetic
