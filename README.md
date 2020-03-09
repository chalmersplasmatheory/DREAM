# TQS
This directory contains the Thermal Quench Simulator (TQS) code.

## Compilation
To compile TQS, go to the root TQS directory and run the following commands:

```bash
$ mkdir -p build
$ cd build
$ cmake ..
$ make -j NTHREADS
```
where ``NTHREADS`` is the number of CPU threads on your computer.

