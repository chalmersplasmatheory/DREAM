Developer information
---------------------
This section covers the design of the DREAM source code. It is intended as an
introduction to the code for new developers of DREAM, with minimal mathematical
details. For an in-depth look at how the equations solved by the code have been
implemented, please consult the DREAM manual.


Conventions
===========
Under this heading we collect information about various non-trivial conventions
used in the code.

Design philosophy
*****************
**Re-usability** (minimize the amount of duplicated code)

Libraries and sub-projects
**************************
The DREAM codebase is divided into three separate components: the FVM library,
the DREAM library, and the DREAMi executable. The purpose is to make the code
as modular and reusable as possible.

The FVM library
^^^^^^^^^^^^^^^
The FVM library implements everything related to the discretisation scheme
(finite volume method) used in DREAM, and can be thought of as a "mathematics"
library (in contrast to the DREAM library, which would be the "physics"
library). The FVM library is completely independent of the other DREAM
components and could thus be re-used in future non-DREAM PDE solvers that
utilize the finite volume method.

The DREAM library
^^^^^^^^^^^^^^^^^
The DREAM library is the kernel of the code. It implements all the physics and
solves all equations of DREAM, and basically everything needed for DREAM
simulations. The only thing the DREAM library doesn't contain is a way for the
user to directly run the code. To run code, one would have to turn to DREAMi,
which is the official DREAM interface, or implement a custom C++ interface.

*Why not make DREAM directly runnable?* The great strength of separating the
kernel code from the user interface is that future, smarter physicists who need
to run DREAM in ways we have not foreseen (e.g. launch thousands of similar
simulations) will easily be able to implement their own interface program which
can talk to the DREAM library while at the same time providing the most
convenient interface for the user. Another, very similar, argument for keeping
the kernel in separate library that we use is that it easily allows us to write
a user-friendly Python interface.

The DREAMi executable
^^^^^^^^^^^^^^^^^^^^^
The DREAMi executable is a very lightweight front-end for DREAM which is only
supposed to pass on information from the user to the library about where DREAM
settings are located. The DREAMi executable is at the time of writing the only
way to run DREAM simulations, but in the future we should be able to provide
other interfaces, most notably a Python interface.

Distribution vs. flux grid
**************************
Since we use a finite volume method for discretising derivatives, two types of
grids appear in all computations: a "distribution" grid, on which the unknown
quantities (such as the distribution function) are determined, and a flux grid,
on which phase space fluxes are evaluated. In the code, one should by default
assume that a quantity is defined on the distribution grid, and if this is not
the case, the variable should be given the suffix ``_f``, indicating that it
is defined on a flux grid. In some cases it will be obvious on which flux
grid the quantity is defined on; there are three to choose among: the radial
flux grid, the first momentum coordinate flux grid, and the second momentum
coordinate flux grid. This is particularly the case for the coordinate vectors,
which are hence called

+-------------------+--------------------------------------------------------+
| **Variable name** | **Description**                                        |
+-------------------+--------------------------------------------------------+
| ``r``             | List of radial coordinates on distribution grid        |
+-------------------+--------------------------------------------------------+
| ``r_f``           | List of radial coordinates on flux grid                |
+-------------------+--------------------------------------------------------+
| ``p1``            | List of first momentum coordinate on distribution grid |
+-------------------+--------------------------------------------------------+
| ``p1_f``          | List of first momentum coordinate on flux grid         |
+-------------------+--------------------------------------------------------+
| ``p2``            | List of second momentum coordinate on distribution grid |
+-------------------+--------------------------------------------------------+
| ``p2_f``          | List of second momentum coordinate on flux grid         |
+-------------------+--------------------------------------------------------+

However, most quantities could be defined on either of the three flux grids
available, and hence we need to indicate which of these it is defined on. A
typical example is the phase space jacobian, denoted by :math:`\mathcal{V}'`.
For the discretisations, we need to know

.. math::

   \texttt{V}_{\tt p} &= \mathcal{V'}\left(r_k, p_{1,i}, p_{2,j}\right),\\
   \texttt{V}_{\tt p\_fr} &= \mathcal{V'}\left( r_{k-1/2}, p_{1,i}, p_{2,j} \right)\\
   \texttt{V}_{\tt p\_f1} &= \mathcal{V'}\left( r_{k}, p_{1,i-1/2}, p_{2,j} \right)\\
   \texttt{V}_{\tt p\_f2} &= \mathcal{V'}\left( r_{k}, p_{1,i}, p_{2,j-1/2} \right)

that is, the phase space jacobian evaluated either on (i) all distribution
grids, (ii) :math:`r` flux grid and :math:`p_1/p_2` distribution grids,
(iii) :math:`p_1` flux grid and :math:`r/p_2` distribution or (iv) :math:`p_2`
flux grid and :math:`r/p_1` distribution grids. To separate between these cases,
add the suffices ``_fr``, ``_f1`` or ``_f2`` to the cases (ii)-(iv) respectively.
One will therefore find the following variables in the code:

+-------------------+--------------------------------------------------------------+
| **Variable name** | **Corresponding to...**                                      |
+-------------------+--------------------------------------------------------------+
| ``Vp``            | :math:`\mathcal{V'}\left(r_k, p_{1,i}, p_{2,j}\right)`       |
+-------------------+--------------------------------------------------------------+
| ``Vp_fr``         | :math:`\mathcal{V'}\left(r_{k-1/2}, p_{1,i}, p_{2,j}\right)` |
+-------------------+--------------------------------------------------------------+
| ``Vp_f1``         | :math:`\mathcal{V'}\left(r_{k}, p_{1,i-1/2}, p_{2,j}\right)` |
+-------------------+--------------------------------------------------------------+
| ``Vp_f2``         | :math:`\mathcal{V'}\left(r_{k}, p_{1,i}, p_{2,j-1/2}\right)` |
+-------------------+--------------------------------------------------------------+


Grid
====
The ``Grid`` object represents a general computational grid, consisting of one
radial coordinate and two momentum coordinates. It is used as an interface for
computing and storing common data related to the grid, such as coordinate
values, grid step lengths and jacobians. Note, though, that actual computations
and storage (with the exception of jacobians) are carried out by the more
specific ``RadialGrid`` and ``MomentumGrid`` classes, which (as their names
suggest) represent the radial and momentum coordinates respectively.

