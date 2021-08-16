# DREAM disruption simulations
This directory contains the script used to configure and run simulations for
section 6 of the DREAM paper (https://doi.org/10.1016/j.cpc.2021.108098). The
simulations are aimed at showcasing the broad range of functionality provided
by DREAM. A detailed walk-through of the example can be found in the DREAM
online documentation: https://ft.nephy.chalmers.se/dream.

## Prerequisites
To run the simulations, you must have installed DREAM on your system. In order
for the ``DREAM.runiface()`` Python function to work you must also have defined
the environment variable ``DREAMPATH`` to point to the DREAM root directory.

## How to run
All configuration (except specific tokamak parameters) and simulation execution
is done in the script ``generate.py``. The script uses the Python ``argparse``
package and, unless you would like to modify simulation settings, can
immediately be run from a terminal with arguments specifying which model and
simulation scenario to run. Running ``./generate.py --help`` yields the
following output:
```
usage: generate.py [-h] [-c] [-e EXTENSION] [--fluid] [--kinetic] [--superthermal]
                   [--isotropic] [--hybrid] [--steady-state-E] [-r RUNFROM] [-s SCENARIO]
                   [-v [VERBOSE ...]]

Run a DREAM disruption simulation

optional arguments:
  -h, --help            show this help message and exit
  -c, --cylindrical     Run in cylindrical geometry
  -e EXTENSION, --extension EXTENSION
                        Append the given string to the end of all file names
  --fluid               Run simulation in pure fluid mode
  --kinetic             Run simulation in fully kinetic mode
  --superthermal        Run simulation in superthermal mode
  --isotropic           Run simulation in isotropic mode
  --hybrid              Run simulation in hybrid superthermal-kinetic mode
  --steady-state-E      Assume that a steady state electric field is applied for initial
                        current density
  -r RUNFROM, --run-from RUNFROM
                        Determines which simulation to start running from (0 = current, 1 =
                        ionization, 2 = CQ
  -s SCENARIO, --scenario SCENARIO
                        Index of scenario to run 0-4
  -v [VERBOSE ...], --verbose [VERBOSE ...]
                        Select which simulations to make verbose (0 = current, 1 = ionization,
                        2 = CQ
```
For example, to generate all data for figures 8 and 9 in the DREAM paper, the
following commands can be run:
```
./generate.py --fluid --scenario 1
./generate.py --isotropic --scenario 1
./generate.py --superthermal --scenario 1
./generate.py --kinetic --scenario 1
```
To generate data for figures 10 and 11, run the following commands:
```
./generate.py --fluid --scenario 4
./generate.py --isotropic --scenario 4
./generate.py --superthermal --scenario 4
./generate.py --kinetic --scenario 4
```
Beware that the superthermal and kinetic simulations can take a significant
amount of time.

