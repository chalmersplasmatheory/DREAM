# DREAM
This directory contains the Disruption Runaway Electron Avoidance Model (DREAM)
code. The **online documentation** is available at https://ft.nephy.chalmers.se/dream.

The official DREAM paper has now been published: [doi:10.1016/j.cpc.2021.108098](https://doi.org/10.1016/j.cpc.2021.108098).
(it is also on arXiv: [2103.16457](https://arxiv.org/abs/2103.16457).)

## Requirements
To compile DREAM, you need to have the following software installed:

- [CMake](https://cmake.org/) >= 3.12
- A C++17 compatible compiler (such as gcc >= 7.0)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/) >= 2.0
- [HDF5](https://www.hdfgroup.org/)
- [PETSc](https://www.mcs.anl.gov/petsc)
- OpenMP
- MPI (for PETSc)
- Python 3 (required for generating ADAS data)

Additionally, to use the DREAM Python interface, you need the following
Python packages:

- h5py
- matplotlib
- numpy
- packaging
- scipy

### Notes on PETSc
While most of the software required by DREAM can be installed directly from
your Linux distribution's package repository, PETSc usually requires a manual
setup. To install PETSc, grab its sources from the PETSc website or clone the
PETSc git repository:
```bash
$ git clone -b release https://gitlab.com/petsc/petsc.git petsc
```
After this, compiling PETSc should be a matter of running the following
commands:
```bash
$ ./configure PETSC_ARCH=linux-c-opt
...
$ make PETSC_DIR=/path/to/petsc PETSC_ARCH=linux-c-opt all
...
```
Optionally, you can also run ``make check`` after ``make all``.

Once PETSc has been compiled with the above commands, you only need to make sure
that DREAM will be able to find your PETSc installation. The easiest way to
achieve this is to add the ``PETSC_DIR`` and ``PETSC_ARCH`` environment
variables used above to your ``~/.bashrc`` file (if you use bash; if you're
unsure, you probably do):
```bash
...
export PETSC_DIR="/path/to/petsc"
export PETSC_ARCH=linux-c-opt
```
Of course, the value for ``PETSC_DIR`` should be modified according to where
you installed PETSc. An alternative to modifying your ``~/.bashrc`` file is to
just give these variables directly to CMake every time you reconfigure DREAM
(which is usually not very often, unless you're a DREAM developer).

## Compilation
To compile DREAM, go to the root DREAM directory and run the following commands:

```bash
$ mkdir -p build
$ cd build
$ cmake ..
$ make -j NTHREADS
```
where ``NTHREADS`` is the number of CPU threads on your computer. If CMake can't
find PETSc, you can change the ``cmake`` command above to read
```bash
$ cmake .. -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH=linux-c-opt
```
where ``/path/to/petsc`` is the path to the directory containing your PETSc
installation.

### "PETSc was configured with MPICH but now appears to be compiling using a non-MPICH mpi.h"
This error can occur if you have only installed MPI locally for PETSc, or if you
have multiple MPI implementations (e.g. MPICH and OpenMPI) installed on your
system. If you installed MPICH automatically during the configuration of PETSc
you should run ``cmake`` with the flag
```bash
$ cmake .. -DMPI_CXX_COMPILER=/path/to/petsc/linux-c-opt/bin/mpicxx
```
Alternatively, if you compiled PETSc with an MPI installation you should specify
```bash
$ cmake .. -DMPI_EXECUTABLE_SUFFIX=.mpich
```
if you compiled PETSc with MPICH, or
```bash
$ cmake .. -DMPI_EXECUTABLE_SUFFIX=.openmpi
```
if you compiled PETSc with OpenMPI.

## Documentation
Online documentation for how to run and extend the code is available at
https://ft.nephy.chalmers.se/dream. LaTeX sources for documentation of the
physics model and various mathematical details can be found under
[doc/notes/](https://github.com/chalmersplasmatheory/DREAM/tree/master/doc/notes).

## Usage in publications
If you use DREAM in your scientific publications, please cite the
[DREAM paper](https://doi.org/10.1016/j.cpc.2021.108098):
```
@article {DREAM,
    title = {DREAM: A fluid-kinetic framework for tokamak disruption runaway electron simulations},
    journal = {Computer Physics Communications},
    volume = {268},
    pages = {108098},
    year = {2021},
    issn = {0010-4655},
    doi = {https://doi.org/10.1016/j.cpc.2021.108098},
    url = {https://www.sciencedirect.com/science/article/pii/S0010465521002101},
    author = {Mathias Hoppe and Ola Embreus and Tünde Fülöp}
}
```

## Authors
DREAM is authored by [Ola Embreus](https://github.com/Embreus)
and [Mathias Hoppe](https://github.com/hoppe93), along with members of the
[Chalmers Plasma Theory group](https://ft.nephy.chalmers.se/). It is currently
maintained by Mathias Hoppe and members of the Chalmers Plasma Theory group.
