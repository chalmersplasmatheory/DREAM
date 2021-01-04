.. _compiling:

Compiling
=========
Building DREAM should be relatively straightforward once all external dependencies
are installed on your system.

Preliminaries
-------------
In order to build DREAM, you must first make sure that the following software is
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

Building DREAM
--------------

Summary
*******
The whole build process is described in more detail below, but can be summarised
using the following chain of commands::

   $ cd /path/to/DREAM
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make -j NTHREADS

where ``NTHREADS`` is the number of CPU threads available on your computer.
(This number if optional, and there is not really any harm in specifying the
"wrong" number; in general, the more "correct" it is, the faster the compilation
will go).


Building PETSc with MUMPS
-------------------------
By far the easiest way of adding MUMPS support to PETSc is by configuring PETSc
with the flags ``--download-mumps`` and ``--download-scalapack`` (ScaLAPACK is
a dependency of MUMPS). Thus, you should configure PETSc in a manner similar to
the following::

    $ ./configure --download-mumps --download-scalapack

If you run into the Fortran error ``Rank mismatch between actual argument at
(1) and actual argument at (2) (scalar and rank-1)``, add the flag
``--FFLAGS=-fallow-argument-mismatch`` to the configure line above.

Building PETSc with Intel MKL
-----------------------------
The Intel Math Kernel Library (MKL) contains the PARDISO linear solver which can
be used by PETSc. To include support for PARDISO in PETSc, configure it with::

    $ ./configure --with-mkl_pardiso-dir=/path/to/mkl --with-blaslapack-dir=/path/to/mkl

where ``/path/to/mkl`` is the path to where the Intel MKL library is installed.
