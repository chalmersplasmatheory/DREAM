.. _runiface:

runiface
========

DREAM simulations can be executed directly from your Python code using the
``runiface()`` function. This function takes either a ``DREAMSettings`` object,
or the name of a file containing such an object, as input and executes the
DREAMi program.

.. contents:: Page overview
   :local:
   :depth: 3

Setting path to DREAM
---------------------
When using ``runiface()``, the DREAM Python interface must be able to locate
the ``dreami`` executable. Normally this is done automatically when the DREAM
Python package is loaded, but if you have an unusual organization of the files
this may fail. In such a case, you can explicitly specify the path to the
``dreami`` executable by setting the ``DREAMPATH`` environment variable.

If you are using the bash shell, add the following line to the file
``~/.bashrc``:

.. code-block:: bash

   export DREAMPATH="/path/to/DREAM"

and make sure to replace ``/path/to/DREAM`` with the actual path to your main
DREAM directory. The DREAM Python interface will then assume that that directory
contains the path ``build/iface``, which in turn should contain the executable
``dreami``.

Using runiface()
----------------
To run DREAM directly from a Python script, simply call ``runiface()`` with a
``DREAMSettings`` object or the name of an HDF5 file containing such an object:

.. code-block:: python

   ds = DREAMSettings()
   ...
   do = runiface(ds)

The ``runiface()`` function will call DREAM and load the resulting output file
into a ``DREAMOutput`` object, which is returned by the function. By default,
DREAM is instructed to store the output in a temporary file which is removed
immediately after being loaded by ``runiface()``. If you want to keep the output
file, you must manually specify the name of the output file to generate:

.. code-block:: python

   # Stores output in 'output.h5', loads the file and returns
   # a DREAMOutput object when finished.
   do = runiface(ds, 'output.h5')

The full syntax of the function is

.. code-block:: python

   do = runiface(ds, 'output.h5', quiet=True, timeout=3600, nthreads=4)

.. note::

   Internally, ``runiface()`` saves the ``DREAMSettings`` to a temporary file
   and runs DREAM with the name of that file as argument. Temporary files with
   weird names may therefore show up while running the code. The files are
   automatically removed after DREAM finishes.

Verbosity
*********
DREAM writes some information to ``stdout`` and ``stderr`` while running. By
default ``runiface()`` forwards this information from DREAM to the terminal. If
you want to suppress this output, you can provide the argument ``quiet=True``
to the function:

.. code-block:: python

   do = runiface(ds, quiet=True)

Timeout
*******
When running many simulations it may be difficult to keep track of whether any
particular simulation has stalled. It is therefore possible to impose a timeout
on the process which causes it to be killed if it runs for longer than a
specified number of **seconds**:

.. code-block:: python

   do = runiface(ds, timeout=True)

Number of threads
*****************
The Intel MKL linear solver allows for parallelization of the matrix inversion
step in each iteration. It has been observed that in certain cases, this may
cause separate DREAM processes to interfere with each other and significantly
slow each other down
(see `Issue #231 <https://github.com/chalmersplasmatheory/DREAM/issues/231>`_).
To avoid this, one can specify the maximum number of threads which each process
is allowed to occupy using the ``nthreads`` argument:

.. code-block:: python

   do = runiface(ds, nthreads=4)

Parallel runiface()
-------------------
A parallel version of ``runiface()`` has also been implemented. This parallel
version is called ``runiface_parallel()`` and is able to launch multiple DREAM
simulations simultaneously. Since the DREAM kernel is not parallelized, this
presents an opportunity to speed up studies requiring multiple independent
simulations to be conducted.

The parallel version of ``runiface()`` takes a list of ``DREAMSettings`` (or
names of files containing such) objects and output filenames as input, both of
which are required arguments. The arguments ``quiet`` and ``timeout`` are also
present in ``runiface_parallel()`` and have the same role on a per-simulation
basis (i.e. if a timeout is imposed, this timeout is enforced for each
individual simulation, and not collectively for all simulations). In addition
to these arguments, it is also possible to specify a list of files to which
``stdout`` and ``stderr`` data will be piped, as well as the number of DREAM
processes which can be kept running at a time, and the number of threads each
process is allowed to use.

Basic usage
***********
The basic use of ``runiface_parallel()`` is illustrated below:

.. code-block:: python

   from DREAM import DREAMSettings

   dss = []
   output = []
   for i in range(10):
       # User-defined function returning a
       # populated DREAMSettings object:
       ds = setup_DREAMSettings(index=i)

       dss.append(ds)
       output.append(f'output{i}.h5')

   # Run simulations in parallel,
   # 4 simulations at a time
   runiface_parallel(dss, output, quiet=True, njobs=4)

Redirecting output
******************
When using ``runiface_parallel()`` all DREAM processes will, by default, print
to the same terminal. The printed information may therefore be very difficult to
interpret, as text from different processes are mixed together. To redirect the
output/error information for each process to separate files, you may use the
``stdout_list`` and ``stderr_list`` arguments to ``runiface_parallel()``:

.. code-block:: python

   from DREAM import DREAMSettings

   dss = []
   output = []
   stdout, stderr = [], []
   for i in range(10):
       # User-defined function returning a
       # populated DREAMSettings object:
       ds = setup_DREAMSettings(index=i)

       dss.append(ds)
       output.append(f'output{i}.h5')
       stdout.append(f'stdout{i}.txt')
       stderr.append(f'stderr{i}.txt')

   # Run simulations in parallel,
   # 4 simulations at a time
   runiface_parallel(dss, output, stdout_list=stdout,
       stderr_list=stderr, quiet=True, njobs=4)

If provided, ``stdout_list`` and ``stderr_list`` must contain equally many
elements as the ``dss`` and ``output`` lists, i.e. one element per DREAM
simulation.

.. note::

   You can follow the progress of each simulation in real time by reading the
   contents of the files. One particularly useful UNIX command for this is
   ``tail``. In the example above, you could follow the progress of the first
   process using the command ``tail -f stdout1.txt``.

Limiting threads
****************
Although the DREAM kernel is not parallelized itself, the linear solver in
Intel MKL (going by the name ``DREAM.Settings.Solver.LINEAR_SOLVER_MKL`` in
DREAM) is able to parallelized using threads. In some circumstances, this
parallelization has be observed to severely impact performance of individual
simulations when using ``runiface_parallel()`` to run multiple simulations in
parallel. The cause of this performance penalty seems to be that the linear
solver uses all available threads on the system, forcing expensive task
switching and causing the processes to run serially in practice. To avoid this,
one can provide the argument ``nthreads`` to ``runiface_parallel()`` with the
maximum number of threads allowed for each process. The value provided to this
argument should be such that ``njobs * nthreads`` equals the number of available
CPU cores, for optimal performance.

.. code-block:: python

   from DREAM import DREAMSettings

   dss = []
   output = []
   for i in range(10):
       # User-defined function returning a
       # populated DREAMSettings object:
       ds = setup_DREAMSettings(index=i)

       dss.append(ds)
       output.append(f'output{i}.h5')

   # Run simulations in parallel,
   # 4 simulations at a time
   runiface_parallel(dss, output, njobs=4, nthreads=2)


