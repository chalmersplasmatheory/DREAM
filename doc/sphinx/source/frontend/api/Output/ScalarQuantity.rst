.. _do-scalarquantity:

ScalarQuantity
==============
The ``ScalarQuantity`` class encapsulates data from scalar quantities in DREAM
output, i.e. data which only has a single, scalar value at a given time. Scalar
quantities only evolve in time. Examples of scalar quantities found in DREAM
are the total plasma current ``I_p`` and the poloidal flux on the tokamak
wall ``psi_edge``.

Usage
-----
.. code-block:: python

   s1 = ScalarQuantity(name='MyScalarQuantity', data=[1,2,3], grid=grid, output=do)
   s2 = ScalarQuantity(name='MyScalarQuantity', data=[4,5,6], grid=grid, output=do)

As with any :ref:`do-unknownquantity`, a ``ScalarQuantity`` can be part of an
arithmetic expression with another ``ScalarQuantity``:

.. code-block:: python

   s3 = s1 + s2
   s4 = s1 - s2
   s5 = s1 * s2
   s6 = s1 / s2

Operations are element-wise.

It is also possible to easily plot a ``ScalarQuantity`` by calling the
``plot()`` member method:

.. code-block:: python

   s1.plot()

If desired, a selection of time points for which to plot the quantity can be
specified to ``plot()``:

.. code-block:: python

   s1.plot(t=[2,3])

Class definition
----------------

.. autoclass:: DREAM.Output.ScalarQuantity.ScalarQuantity
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __getitem__

