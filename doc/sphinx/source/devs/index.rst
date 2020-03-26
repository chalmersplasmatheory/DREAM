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
| ``p2``            | List of first momentum coordinate on distribution grid |
+-------------------+--------------------------------------------------------+
| ``p2_f``          | List of first momentum coordinate on flux grid         |
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
values, grid step lengths and jacobians. Note, though, that actual copmutations
and storage (with the exception of jacobians) are carried out by the more
specific ``RadialGrid`` and ``MomentumGrid`` classes, which (as their names
suggest) represent the radial and momentum coordinates respectively.

