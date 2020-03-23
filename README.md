# DREAM
This directory contains the Disruption and Runaway Electron Avoidance Model (DREAM) code.

## Compilation
To compile DREAM, go to the root DREAM directory and run the following commands:

```bash
$ mkdir -p build
$ cd build
$ cmake ..
$ make -j NTHREADS
```
where ``NTHREADS`` is the number of CPU threads on your computer.

