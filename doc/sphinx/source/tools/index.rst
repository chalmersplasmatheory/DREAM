.. _dream-tools:

DREAM tools
===========
DREAM comes with a number of helpful tools which can be useful when working
with DREAM input data.

.. toctree::
   :maxdepth: 2

   atomicdata
   eqget

Atomic data
-----------
Adding support for new atomic elements is relatively straightforward, but
requires some manual work. On the page about Atomic

Read about :ref:`dream-atomicdata`.


Debugging tools
---------------
.. todo::

   Tools for advanced debugging of DREAM simulations.


Equilibrium tools
-----------------
DREAM provides some Python classes and graphical tools for working with magnetic
equilibria in various formats. As described on the page about the
:ref:`_radialgrid-numerical`, DREAM uses a LUKE format for equilibrium data.
We however also provide a number of classes for converting data from other
formats to the LUKE format.

Read about :ref:`dream-eqget`.


CORSICA profiles reader
-----------------------
.. todo::

   Describe the script ``tools/CORSICAProfileReader.py``.


Transport coefficient reader
----------------------------
.. todo::

   Describe the script ``tools/TransportCoefficientReader.py``.


