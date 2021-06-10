.. _python-frontend:

Python frontend
===============
This section of the documentation deals with the Python frontend built for
running and interacting with DREAM and its data. The code comprising the Python
frontend is located under ``py/`` and comprises three separate projects:

- ``py/DREAM/``: the DREAM Python API which contains all scripts used to
  programmatically interact with DREAM input and output.
- ``py/cli/cli.py``: simple script for interactively working with DREAM output.
- ``py/theater``: The DREAM GUI interface which can be used to both configure
  and run DREAM simulations. **Not yet implemented**

Python API
----------
The Python API, located under ``py/DREAM/``, is the engine that powers all parts
of the Python interface. It does not contain any runnable scripts, but rather
holds several helper classes and functions which significantly simplify
interaction with DREAM input and output. Aside from a few helper functions and
classes located directly under the ``DREAM`` namespace, the interface is
divided into two sections: ``DREAM.Settings`` and ``DREAM.Output``. The former
contains all code concerning DREAM input, while the latter deals with the
output.

.. toctree::
   :maxdepth: 2

   api/ConvergenceScan
   api/ConvergenceScanPlot
   api/DREAMOutput
   api/Settings/DREAMSettings
   api/runiface

Interactive CLI
---------------
The script ``py/cli/cli.py`` facilitates easier access to and interaction with
DREAM output files. The script runs in interactive mode (i.e. ``python3 -i``)
and loads either the specified file or the file ``output.h5`` in the current
working directory as a :ref:`DREAMOutput` object. Additionally, all unknowns of
the equation system are declared as global variables so that they can be
directly access by name in the Python REPL. Example session:

.. code-block:: bash

   $ /path/to/DREAM/py/cli/cli.py my-dream-output.h5
   Loaded 8 unknowns (9.21 MiB)
   DREAM Grid Object
      t: 11 elements from 0.0 to 0.010000000000000002
      r: 10 elements from 0.011 to 0.20899999999999996
      hottail   (p/xi)
         p:     1000 elements from 0.0005 to 0.9995
         xi:    10 elements from -0.9 to 0.9
      runaway   (DISABLED)

   Unknowns:
      E_field, f_hot, n_cold, n_hot, n_i, n_re, n_tot, T_cold
   >>> f_hot
   (f_hot) Kinetic quantity of size NT x NR x NXI x NP = 11 x 10 x 10 x 1000
   >>> np.sum(f_hot[:])
   3.493714926993402e+26
   >>> exit()

After running through the ``cli.py`` script which initializes the session, the
Python interpreter switches to interactive mode allows the user to interact as
usual. Everything you can normally do in the Python REPL can be done after
running ``cli.py``.

Most of the logic of the ``cli.py`` script is intentionaly kept in the function
``setup_interactive`` of the Python API. This allows users to easily put
together their own versions of ``cli.py`` specialized to their particular
needs. For example, the following variation of ``cli.py`` sets up an interactive
environment and imports all ``numpy`` functions into the global namespace:

.. code-block:: python

   import matplotlib.pyplot as plt
   from numpy import *

   from DREAM import setup_interactive
   from DREAM.DREAMOutput import DREAMOutput

   do = DREAMOutput('output.h5')
   setup_interactive(do, glob=globals())

   # From this point on, the user will be able to interact
   # with the DREAM output in the Python REPL...

All variables available in the interactive session are of the types defined in
the Python API. More specifically, the unknowns of a DREAM Output equation
system are loaded with the following types (at the time of writing):

+---------------+----------------------------------------+
| Variable name | Data type                              |
+===============+========================================+
| ``E_field``   | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``f_hot``     | :ref:`do-hotelectrons`                 |
+---------------+----------------------------------------+
| ``f_re``      | :ref:`do-reelectrons`                  |
+---------------+----------------------------------------+
| ``I_p``       | :ref:`do-scalarquantity`               |
+---------------+----------------------------------------+
| ``j_hot``     | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``j_ohm``     | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``j_tot``     | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``n_cold``    | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``n_hot``     | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``n_i``       | :ref:`do-ionhandler`                   |
+---------------+----------------------------------------+
| ``n_re``      | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``n_tot``     | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``psi_edge``  | :ref:`do-scalarquantity`               |
+---------------+----------------------------------------+
| ``psi_p``     | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``T_cold``    | :ref:`do-fluidquantity`                |
+---------------+----------------------------------------+
| ``x_p``       | :ref:`do-spishardpositions`            |
+---------------+----------------------------------------+
| ``Y_p``       | :ref:`do-spishardradii`                |
+---------------+----------------------------------------+

All other variables are loaded as :ref:`do-unknownquantity` if the type has not
been explicitly specified in the dict ``SPECIAL_TREATMENT`` in
:ref:`do-eqsys`.

The DREAM Theater (GUI)
-----------------------
*The DREAM Theater GUI has not yet been implemented*
