ConvergenceScan
===============
The ``ConvergenceScan`` class provides a simple way of quickly generating
parameter scans, particularly in the resolution parameters. A
``ConvergenceScan`` object holds a reference to a baseline :ref:`DREAMSettings`
object, as well as a list of input and output parameters. Each input and output
parameter has a name associate with it, as well as a function which modifies a
given :ref:`DREAMSettings` object (input) or extracts data from a
:ref:`DREAMOutput` object (output).

Once the ``ConvergenceScan`` has been successfully executed, the result can be
conveniently plotted using the :ref:`ConvergenceScanPlot` class.

The most basic use case is illustrated by the following example (although far
more complicated things can be achieved with the class):

.. code-block:: python

   import DREAM
   # ... or
   # from DREAM.ConvergenceScan import ConvergenceScan
   # ... and use as just "ConvergenceScan()"

   # First, create a DREAMSettings object as you normally would...
   ds = DREAM.DREAMSettings()

   # set up ds...
   ...

   # Create convergence scan
   cs = DREAM.ConvergenceScan(settings=ds,
        inparams=['nt', 'hottail.np', 'hottail.nxi'],
        outparams=['other.fluid.runawayRate'])

   # Run convergence scan 
   cs.run()

   # Save results to file
   cs.save('convergence.h5')


.. note:: **DREAM path**

   In order for the convergence scan to find DREAM, you must first properly set
   up your environment before using the ``ConvergenceScan`` object. Primarily,
   the environment variable ``DREAMPATH`` variable should be defined and point
   to the DREAM source directory. You can read more about how to do this on the
   page :ref:`runiface`.


.. _convergencescan-by-name:

Input/output parameter by name
------------------------------
Input and output parameters can be specified by name. Except for a few special
parameters, the name is the full name, including parent class names in the
owning :ref:`DREAMSettings`/:ref:`DREAMOutput` object. For example, to vary the
radial grid resolution, which is accessed in the :ref:`DREAMSettings` object
using ``DREAMSettings.radialgrid.nr``, one would generally specify the name of
the parameter as ``radialgrid.nr``, that is, everything after the the
``DREAMSettings.``. The same applies to :ref:`DREAMOutput`, so that the name
``eqsys.I_p`` would access the final plasma current.

Parameters can be added either when the ``ConvergenceScan`` object is created,
or afterwards (or both):

.. code-block:: python

   ds = DREAMSettings()
   ...

   cs = ConvergenceScan(ds, inparams=['nt'], outparams=['other.fluid.runawayRate'])
   cs.addScanParameter('radialgrid.nr')
   cs.addScanParameter('hottailgrid.pgrid.np')
   cs.addScanParameter('hottailgrid.pgrid.nxi')

Aliases for input parameters
****************************
Some input parameters are special and have aliases which allow for easier
access. These are listed in the table below:

+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| Full name                                                | Alias                 | Description                                                                |
+==========================================================+=======================+============================================================================+
| ``hottailgrid.pgrid.np``                                 | ``hottailgrid.np``    | Momentum resolution on hottail grid.                                       |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``hottailgrid.xigrid.nxi``                               | ``hottailgrid.nxi``   | Pitch resolution on hottail grid.                                          |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``hottailgrid.pgrid.np`` and ``hottailgrid.xigrid.nxi``  | ``hottail``           | Momentum AND pitch resolution on hottail grid (adds two input parameters). |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``radialgrid.nr``                                        | ``nr``                | Radial resolution.                                                         |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``timestep.nt``                                          | ``nt``                | Time resolution.                                                           |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``runawaygrid.pgrid.np``                                 | ``runawaygrid.np``    | Momentum resolution on runaway grid.                                       |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``runawaygrid.xigrid.nxi``                               | ``runawaygrid.nxi``   | Pitch resolution on runaway grid.                                          |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+
| ``runawaygrid.pgrid.np`` and ``runawaygrid.xigrid.nxi``  | ``runaway``           | Momentum AND pitch resolution on runaway grid.                             |
+----------------------------------------------------------+-----------------------+----------------------------------------------------------------------------+

Special notes on output parameters
**********************************
Output parameter data is accessed by index, meaning that any
:ref:`do-unknownquantity` or :ref:`do-otherquantity` may be used as the output
parameter. When the output parameter is given by name, the code will access only
the very last element of the data. This means that, for

- :ref:`do-scalarquantity`'s (such as ``I_p``), the final value is used.
- :ref:`do-fluidquantity`'s (such as ``E_field``), the final value at the outermost radius is used.
- :ref:`do-kineticquantity`'s (such as ``f_hot``), the final value at the outermost radius in :math:`\xi=1` and :math:`p=p_{\rm max}` is used.

When specifying output parameters by name, the very last element of the 
parameter data will generally be accessed.

Input parameter custom function
-------------------------------
Sometimes the parameter to scan is more complicated to set than just
increasing/decreasing it by a constant float value. In this case, one can
instead define a custom function which modifies the settings object in a
successive fashion. The function can only be provided via a call to
``addScanParameter()`` and must thus be added after the ``ConvergenceScan``
object has been constructed.

The ``addScanParameter()`` should be called in the following way:

.. code-block:: python

   cs.addScanParameter(name='inparam', f=customFunction, baselineValue=baseval)

These are the required parameters; the other parameters can also be specified if
desired.

The name of the parameter is essentially arbitrary. It is only used by the
``ConvergenceScan`` object to identify the baseline value if the
``baselineValue`` parameter is **not** specified. If this is the case, the
baseline value is taken from the variable in the previously given
:ref:`DREAMSettings` object which has the given name.

The custom function ``customFunction()`` should have the following signature:

.. py:function:: customInputFunction(index, settings, baseline)

   :param int index:              Index of simulation to set up (``0`` means the baseline case; negative values are possible)
   :param DREAMSettings settings: Settings object to modify. On input, this object is a copy of the baseline settings object specified when constructing the ``ConvergenceScan`` object.
   :param baseline:               Baseline value for this parameter.
   :return:                       Tuple consisting of the modified settings object and a ``float`` representing the value set.

.. note:: **Lambda expressions**

   The use of lambda functions is often appropriate when passing functions
   to ``addScanParameter()`` and can provide more compact code
   (see, for example,
   https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions).

The function is supposed to modify the parameter ``settings``, which is a copy
of the :ref:`DREAMSettings` object given to the ``ConvergenceScan`` when
constructed, and return the modified settings object along with a numerical
value representing the assigned setting (even if the value is not numerical
itself; then it could for example be ``index``). An example implementation is:

.. code-block:: python

   def _CS_getiNt(index: int, settings: DREAMSettings, baseline):
       # Calculate new value to set
       val = max(1,int(np.round(baseline * np.float_power(2, index))))
       # Modify settings object
       settings.timestep.setNt(val)
       # Return modified object and new parameter value
       return settings.val

Simulation indices
******************
The ``index`` parameter accepted by the custom function indicates the stage of
the convergence scan to set up. Indices work such that ``0`` correspond to the
baseline case, while positive values indicate "higher resolution" and negative
values indicate "lower resolution" (of course, users are welcome to
reinterpret the distinction between positive and negative indices however they
desire). The scanner expects the function to modify the object in a
deterministic way so that a call with a specific index always results in the
same settings being applied. In general, the baseline case is only run once,
instead of once for each scan parameter.

By default, the starting index is ``-1``, which is then gradually increased
until the upper index limit is reach, which is set to ``1`` (inclusive) by
default.

Output parameter custom function
--------------------------------
As with complicated input parameters, more complicated output parameters can
also be accessed via a custom function. The custom function has to be added
separately via the ``addOutputParameter()`` function in the following way:

.. code-block:: python

   cs.addOutputParameter(name='outparam', f=customFunction)

Optionally, a relative tolerance used for the continuous convergence scan mode
can also be provided, but is not required.

If a custom function is provided, the name of the output parameter is only used
when communicating with the user and has no internal significance.

The custom function should have the signature

.. py:function:: customOutputFunction(do: DREAMOutput)

   :param DREAMOutput do: Output object to extract parameter value from.
   :return:               The value of the output parameter.
   :rtype:                float

.. note:: **Lambda expressions**

   The use of lambda functions is often appropriate when passing functions
   to ``addOutputParameter()`` and can provide more compact code
   (see, for example,
   https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions).

The purpose of the function is to process the given :ref:`DREAMOutput` object
in order to obtain the value of the output parameter resulting from the
simulation. The returned value must be a ``float``, although the actual value
returned generally is of little or no interest to the ``ConvergenceScan``
object (it could for example be a binary value, varying discretely between
0 and 1). The only time the ``ConvergenceScan`` object makes a decision based
on the value of the parameter is when ``scanUntilConvergence`` is set to
``True`` for an input parameter.

An example implementation of the custom output parameter function is the
following:

.. code-block:: python

   def customFunction(do: DREAMOutput) -> float:
       # Calculate kinetic energy carried by hot electrons
       Wk = do.eqsys.f_hot.kineticEnergy()[-1,:]
       # Turn into a FluidQuantity
       Wk = DREAM.Output.FluidQuantity('Wk', Wk, do.grid, do)

       # Return total final kinetic energy
       return Wk.integral(t=-1)

Scanning until convergence
--------------------------
If you would like to find the point of convergence without manually tweaking
the scan, you can set the parameter ``scanUntilConvergence=True`` when calling
``addScanParameter()``. This will cause the ``ConvergenceScan`` to increase the
simulation index by one until the output parameters vary by less than the
relative tolerance specified for each output parameter.

.. warning::

   Note that there is a maximum number of permitted iterations in the
   ``ConvergenceScan`` object, called ``NMAX``. It defaults to ``10``, meaning
   that if convergence has not been reached on the tenth iteration
   (corresponding to index ``8`` with the default starting index ``-1``), then
   the scan will stop with an error.

Class documentation
-------------------

.. autoclass:: DREAM.ConvergenceScan
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

