.. _cpp-kernel:

The C++ kernel
==============
In this section of the documentation you will find information about the
computational kernel of DREAM. The kernel consists of two C++ libraries, named
``libfvm`` and ``libdream``, which implement all physics and key algorithms of
the code. While these libraries constitute most of what we consider to be
"the code", the average user cannot directly interact with these libraries but
must rather do so through either the ``dreami`` interface, and/or the Python
interface. This documentation is therefore intended for people who will work
directly with code in the computational kernel of DREAM, and *not* for people
who merely intend to use the code for simulations.

**Separation into two libraries:**
In order to keep the code modular and allow some code to be reused in future,
unrelated applications, the kernel code is separated into two libraries, called
``libdream`` and ``libfvm``. The distinction made between what goes into which
library is somewhat arbitrary, but the general idea is that code which is very
specific to the operation of the DREAM code specifically should go into
``libdream``, while ``libfvm`` should contain code which could be used by any
PDE solver. This distinction can also be phrased as *"libdream implements
the physics, while libfvm implements the mathematics"*. Perhaps more
importantly, ``libfvm`` should be completely independent of any other DREAM
module, whereas ``libdream`` necessarily depends heavily on ``libfvm``.

DREAM --- the physics library
-----------------------------
The DREAM library implements the actual physics simulated by DREAM. It also
contains the C++ API necessary to run simulations. If you would like to build
your own interface to DREAM, then you should interface with ``libdream``.

In addition to implementing calculations of collision frequencies, ionization
rates and the Fokker-Planck equation, ``libdream`` contains the linear and
non-linear solvers, the time stepping algorithms, as well as the initialization
logic.

**Contents**

.. toctree::
   :maxdepth: 3

   dream/index.rst


FVM --- the mathematics library
-------------------------------
The FVM library implements various helper classes which are used by the DREAM
library. Perhaps most importantly, the FVM library implements the ``Grid`` class
and its children, as well as discretizations for the finite volume method (which
is the source of the library's name).

While ``libfvm`` is an essential dependency for the DREAM library, ``libfvm``
itself is intended to note depend on any other DREAM modules. This allows its
code to be reused in future, non-DREAM related applications.

**Contents**

.. toctree::
   :maxdepth: 2

   fvm/index.rst

