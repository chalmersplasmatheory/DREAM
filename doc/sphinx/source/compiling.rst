.. _compiling:

Compiling
=========
Building DREAM should be relatively straightforward once all external dependencies
are installed on your system.

.. contents:: Page overview
   :local:
   :depth: 3

Preliminaries
-------------
In order to build DREAM, you must first make sure that the following software is
installed on your system:

- `CMake <https://cmake.org/>`_ version 3.12 or later
- A C++17 compliant compiler and standard library (i.e. `gcc` >= 7 or equivalent)
- `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_ version 2.0 or later
- `HDF5 <https://www.hdfgroup.org>`_
- `PETSc <https://www.mcs.anl.gov/petsc/index.html>`_
- OpenMP
- MPI (for PETSc)
- Python 3 (required for generating ADAS data)

If you are running Linux, these tools are usually available directly in the
package repository for your Linux distribution. On distributed computers, you
can usually install these programs using a module system (e.g.
`module <http://modules.sourceforge.net/>`_).

Additionally, the SOFT support library
`softlib <https://github.com/hopp93/softlib>`_ is required, but is automatically
downloaded and installed if you have `git <https://git-scm.com/>`_ installed.

Python
******
To use the DREAM Python interface, you will need Python 3 installed on your
system. Additionally, the following packages are needed (all of which should
be installable via ``pip`` (https://pypi.org/):

- h5py
- matplotlib
- numpy
- packaging
- scipy

**Installing dependencies with pip:**

.. code-block:: bash

   $ pip install h5py matplotlib numpy packaging scipy

Building PETSc
--------------
.. tip::

   If you already have an installation of PETSc available (as is the case on
   many public computer systems) it should not be necessary to rebuild PETSc.
   Instead, you can proceed to building DREAM directly.

DREAM uses PETSc, the *Portable, Extensible Toolkit for Scientific Computation*,
for working with sparse matrices. Much of the execution time is spent in PETSc
or the external packages it calls during a DREAM simulation, and it is therefore
useful to pay extra attention to the configuration of the library.

General steps
*************
PETSc can be downloaded using ``git``:

.. code-block:: bash

   $ git clone -b release https://gitlab.com/petsc/petsc.git petsc

We generally recommend using the latest version available of PETSc (which the
above command will give you). The last argument to ``git`` above, ``petsc``,
specifies the location to download PETSc to. As specified above, PETSc is
downloaded to a subdirectory called ``petsc`` located in the current working
directory. A common place to put PETSc is either in your home directory
(``~/petsc``) or, if you are the system administrator, in ``/opt/petsc``.

Once PETSc is downloaded, you will need to configure PETSc. This is achieved
by going into the PETSc directory and running ``./configure``:

.. code-block:: bash

   $ ./configure

.. note::

   Generally you will also want to specify a number of configuration flags.
   See for example :ref:`compiling-recommended-petsc` and
   :ref:`compiling-petsc-external` for details.

After configuration has finished successfully, you can compile using the command

.. code-block:: bash

   $ make all

After this command finishes successfully, you can proceed to compiling DREAM.

.. tip::

   To reduce the amount of typing when compiling DREAM, you can export
   appropriate values for the ``PETSC_DIR`` and ``PETSC_ARCH`` environment
   variables in your ``~/.bashrc`` file:

   .. code-block:: bash

      ...
      export PETSC_DIR="/path/to/petsc"
      export PETSC_ARCH=linux-c-opt

   The values to use for ``PETSC_DIR`` and ``PETSC_ARCH`` are given at the end
   of the PETSc configuration.

.. _compiling-recommended-petsc:

Recommended configuration
*************************
We recommend configuring PETSc with the following command (assuming a GCC or
compatible compiler is used to compile PETSc):

.. code-block:: bash

   $ ./configure --with-debugging=0 --COPTFLAGS="-O3 -march=native -mtune=native" --CXXOPTFLAGS="-O3 -march=native -mtune=native" --FOPTFLAGS="-O3 -march=native -mtune=native"

Of course, you may need additional flags specific for your system, and if you
want support for external solver packages (which we highly recommend!) you will
also need to append the flags described below.

In the suggested line above, the ``--with-debugging=0`` flag disables all debug
settings in PETSc and allows compilation with optimizations. The
``--COPTFLAGS``, ``--CXXOPTFLAGS`` and ``--FOPTFLAGS`` specifies additional
optimization flags to be passed on to the C, C++ and Fortran compilers while
building PETSc. In this case, we use the highest optimization level (``-O3``)
and allow the compiler to use CPU instructions specifically for the system that
PETSc is being compiled on (``-march=native`` and ``-mtune=native``).

.. warning::

   If you are compiling PETSc on an architecturally different system than you
   are going to run the code on, you should not use the ``-march=native`` and
   ``-mtune=native`` flags as this may cause the compiler to generate invalid
   code for the system on which the program will run.

.. _compiling-petsc-external:

External packages
*****************
PETSc provides an interface a large number of external linear solver packages.
In DREAM we have added explicit support for a few of them, and we generally
recommend using one of the packages over the default built-in PETSc sparse LU
factorization algorithm. In our experience, the fastest and most reliable linear
solver when used in conjunction with DREAM is **Intel MKL's PARDISO** solver.

Intel MKL PARDISO
^^^^^^^^^^^^^^^^^
The Intel Math Kernel Library (MKL) contains the PARDISO linear solver which can
be used by PETSc. To include support for PARDISO in PETSc, configure it with::

    $ ./configure --with-mkl_pardiso-dir=/path/to/mkl --with-blaslapack-dir=/path/to/mkl

where ``/path/to/mkl`` is the path to where the Intel MKL library is installed.

The solver can be installed along with the rest of Intel MKL and is available
in the package repositories of many popular Linux distributions (including
Ubuntu 20.04+ and Arch Linux). To install on recent versions of Ubuntu, simply
run

.. code-block:: bash

   sudo apt install intel-mkl

If Intel MKL is *not* available in the package repositories of your Linux
distribution, you can download it from the official
`Intel MKL website <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html>`_.

SuperLU
^^^^^^^
To add support for the SuperLU linear solver to PETSc, configure with the
command::

    $ ./configure --download-superlu


MUMPS (not recommended)
^^^^^^^^^^^^^^^^^^^^^^^
.. warning::

   Although DREAM has support for running with MUMPS, we have experienced
   several stability issues with the MUMPS solver. It is among the fastest
   solvers available for DREAM, but can sometimes fail to invert the equation
   system.

By far the easiest way of adding MUMPS support to PETSc is by configuring PETSc
with the flags ``--download-mumps`` and ``--download-scalapack`` (ScaLAPACK is
a dependency of MUMPS). Thus, you should configure PETSc in a manner similar to
the following::

    $ ./configure --download-mumps --download-scalapack

If you run into the Fortran error ``Rank mismatch between actual argument at
(1) and actual argument at (2) (scalar and rank-1)``, add the flag
``--FFLAGS=-fallow-argument-mismatch`` to the configure line above.

Troubleshooting
***************

Building DREAM
--------------

Summary
*******
The whole build process is described in more detail below, but can be summarised
using the following chain of commands:

.. code-block:: bash

   $ cd /path/to/DREAM
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make -j NTHREADS

where ``NTHREADS`` is the number of CPU threads available on your computer.
(This number if optional, and there is not really any harm in specifying the
"wrong" number; in general, the more "correct" it is, the faster the compilation
will go).

If the ``PETSC_DIR`` and ``PETSC_ARCH`` environment variables are not exported
in your ``~/.bashrc``, you may also need to give them explicitly to CMake:

.. code-block:: bash

   $ cmake .. -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=linux-c-opt

(note that in order for the variables to be defined, you must restart ``bash``
after exporting them in your ``~/.bashrc``).

Troubleshooting
***************

PETSC_EXECUTABLE_RUNS missing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
On some systems, particularly Ubuntu, you will need to override the
``PETSC_EXECUTABLE_RUNS`` CMake variable:

.. code-block:: bash

   $ cmake .. -DPETSC_EXECUTABLE_RUNS=YES

"PETSc was configured with MPICH but now appears to be compiling using a non-MPICH mpi.h"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This error can occur if you have installed MPI while configuring PETSc, or if
you have multiple MPI implementations (e.g. MPICH and OpenMPI) installed
alongside each other on your system. If you installed MPICH automatically during
the configuration of PETSc you should run CMake with the flag

.. code-block:: bash

   $ cmake .. -DMPI_CXX_COMPILER=/path/to/petsc/$PETSC_ARCH/bin/mpicxx

Alternatively, if you compiled PETSc with a system-wide MPICH installation you
should specify

.. code-block:: bash

   $ cmake .. -DMPI_EXECUTABLE_SUFFIX=.mpich

or, you use OpenMPI

.. code-block:: bash

   cmake .. -DMPI_EXECUTABLE_SUFFIX=.openmpi

