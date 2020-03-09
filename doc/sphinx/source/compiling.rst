.. _compiling:

Compiling
=========
Building TQS should be relatively straightforward once all external dependencies
are installed on your system.

Preliminaries
-------------
In order to build TQS, you must first make sure that the following software is
installed on your system:

- `CMake <https://cmake.org/>`_ version 3.9 or later
- A C++17 compliant compiler and standard library (i.e. `gcc` >= 7 or equivalent)
- `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_ version 2.0 or later
- `HDF5 <https://www.hdfgroup.org>`_
- `PETSc <https://www.mcs.anl.gov/petsc/index.html>`_

If you are running Linux, these tools are usually available directly in the
package repository for your Linux distribution. On distributed computers, you
can usually install these programs using a module system (e.g.
`module <http://modules.sourceforge.net/>`_

Additionally, the SOFT support library
`softlib <https://github.com/hopp93/softlib>`_ is required, but is automatically
downloaded and installed if you have `git <https://git-scm.com/>`_ installed.

Building TQS
------------

Summary
*******
The whole build process is described in more detail below, but can be summarised
using the following chain of commands::

   $ cd /path/to/TQS
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make -j NTHREADS

where ``NTHREADS`` is the number of CPU threads available on your computer.
(This number if optional, and there is not really any harm in specifying the
"wrong" number; in general, the more "correct" it is, the faster the compilation
will go).


