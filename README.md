# DREAM

![The DREAM logo](https://github.com/chalmersplasmatheory/media/logo1.png "The DREAM logo")

This directory contains the Disruption Runaway Electron Analysis Model (DREAM)
code. The **online documentation** is available at https://ft.nephy.chalmers.se/dream.

DREAM is a physics simulation framework developed for studying relativistic
runaway electrons in [tokamak](https://en.wikipedia.org/wiki/Tokamak) fusion
devices. Specifically, DREAM solves a system of non-linear partial differential
equations which describe the time evolution of a tokamak plasma. What sets DREAM
apart from other tokamak transport codes is its wide range of models for
studying runaway electron generation and dynamics. In particular, the fluid
equations solved by DREAM can be coupled to a set of kinetic equations for
electrons in order to more accurately describe the runaway electrons.

The official DREAM paper is
[doi:10.1016/j.cpc.2021.108098](https://doi.org/10.1016/j.cpc.2021.108098)
(it is also on arXiv: [2103.16457](https://arxiv.org/abs/2103.16457)).

## Requirements
To compile DREAM, you need to have the following software installed:

- [CMake](https://cmake.org/) >= 3.12
- A C++17 compatible compiler (such as gcc >= 7.0)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/) >= 2.4
- [HDF5](https://www.hdfgroup.org/)
- [PETSc](https://www.mcs.anl.gov/petsc)
- Python 3 (required for generating ADAS data and using the Python interface)

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
$ ./configure PETSC_ARCH=linux-c-opt --with-mpi=0
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
The value for ``PETSC_DIR`` should be modified according to where you installed
PETSc. An alternative to modifying your ``~/.bashrc`` file is to just give these
variables directly to CMake every time you reconfigure DREAM (which is usually
not very often, unless you're a DREAM developer).

## Compilation
*If you're trying to install DREAM on a cluster which has previously been used
to run DREAM, have a look in the
[setup](https://github.com/chalmersplasmatheory/DREAM/tree/master/setup)
directory to see if a build script is already available.*

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

## Documentation
Online documentation for how to run and extend the code is available at
https://ft.nephy.chalmers.se/dream. LaTeX sources for documentation of the
physics model and various mathematical details can be found under
[doc/notes/](https://github.com/chalmersplasmatheory/DREAM/tree/master/doc/notes).

## Citing DREAM
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
    doi = {10.1016/j.cpc.2021.108098},
    url = {https://doi.org/10.1016/j.cpc.2021.108098},
    author = {Mathias Hoppe and Ola Embreus and Tünde Fülöp}
}
```
If you use certain functionality, you should also cite the relevant reference
for that functionality:

- Shattered Pellet Injection: [O. Vallhagen, MSc thesis, Chalmers University of Technology (2021)](https://hdl.handle.net/20.500.12380/302296)
- Plasmoid drift of SPI shards: [O. Vallhagen *et al*, JPP **89** 905890306 (2023)](https://doi.org/10.1017/S0022377823000466)
- Fluid runaway ionization: [M. Hoppe *et al*, accepted for publication in PPCF (2025)](https://arxiv.org/abs/2412.14721)
- Runaway scrape-off: [O. Vallhagen *et al*, accepted for publication in JPP (2025)](https://arxiv.org/abs/2410.03512)

## Development
DREAM development is overseen by the DREAM Developer Council which consists of

![The DREAM Developer Council 2024](https://github.com/chalmersplasmatheory/media/DDC/Council-2024.jgp "The DREAM Developer Council in 2024)

- [Mathias Hoppe](https://www.kth.se/profile/mhop?l=en), KTH Royal Institute of Technology, Stockholm, Sweden
- [Ida Ekmark](https://ft.nephy.chalmers.se/?p=people&id=61), Chalmers University of Technology, Gothenburg, Sweden
- [Lorenzo Votta](https://www.kth.se/profile/votta?l=en), KTH Royal Institute of Technology, Stockholm, Sweden
- [Oskar Vallhagen](https://ft.nephy.chalmers.se/?p=people&id=16), Chalmers University of Technology, Gothenburg, Sweden
- [Peter Halldestam](https://www.ipp.mpg.de/person/140267/5497846), Max Planck Institute for Plasma Physics, Garching-bei-München, Germany

The code was originally developed within the
[Plasma Theory group at Chalmers](https://ft.nephy.chalmers.se/) but is now
coordinated from KTH Royal Institute of Technology, still in close collaboration
with Chalmers. Over the lifetime of the code, a large number of people from
across the world have contributed to its development. A regularly updated list
of contributors can be found in
[CONTRIBUTORS.md](https://github.com/chalmersplasmatheory/DREAM/blob/master/CONTRIBUTORS.md).

The archive of DREAM newsletters, where changes to DREAM are recorded, can be
found at [https://ft.nephy.chalmers.se/dreamnews/](https://ft.nephy.chalmers.se/dreamnews/).

The original idea for DREAM was proposed by
[Ola Embreus](https://github.com/Embreus) who, together with Mathias Hoppe,
developed the first version.
