.. _ds-momentumgrid:

MomentumGrid
=============
The ``hottailgrid`` and ``runawaygrid`` sections of a ``DREAMSettings`` object
specify settings for the hot-tail and runaway momentum grids. Both grids are of
the type ``MomentumGrid``, have the same options and grid types. The only
difference between the grids appears in the kernel where the runaway grid has
:math:`p^{\rm RE}_{1/2}=0` if the hot-tail grid is *disabled*, and
:math:`p^{\rm RE}_{1/2}=p^{\rm HT}_{\rm max}` if *enabled* (where
:math:`p^{\rm HT}_{\rm max}` is the maximum momentum on the hot-tail grid).

.. note::

   Since DREAM uses the finite-volume method, each momentum grid consists of
   two shifted momentum grids referred to as the *flux* and *distribution*
   grids.

   The **flux grid** is indexed with half-indices, starting with :math:`1/2`
   which denotes the point on the lower boundary, and ending with :math:`N+1/2`
   or :math:`\rm max`, denoting the point on the upper boundary. The flux grid 
   is so called because it is where particle fluxes are evaluated.

   The points of the **distribution grid** are situated in between the points
   of the flux grid. They denote the points in which the distribution function
   is calculated.

Momentum coordinates
--------------------
As of this writing, DREAM can only use spherical coordinates (i.e. momentum
magnitude :math:`p` and cosine of pitch angle :math:`\xi`). The code has however
been designed with the intention of allowing the use of cylindrical coordinates
at some point in the future (i.e. momentum parallel :math:`p_\parallel` and
perpendicular :math:`p_\perp` to the magnetic field).

Object documentation
--------------------
.. py:class:: MomentumGrid

.. py:method:: MomentumGrid(name, enabled=True, ttype=MOMENTUMGRID_TYPE_PXI, np=100, nxi=1, pmax=None)

   Construct a new MomentumGrid object.

   :param str name: Name of the momentum grid (either ``hottailgrid`` or ``runawaygrid``
   :param bool enabled: Whether or not the grid is to be used during the simulation.
   :param int type: Type of momentum grid (currently, only ``MOMENTUMGRID_TYPE_PXI``, for spherical coordinates, is supported).
   :param int np: Number of distribution grid points in the spherical coordinate :math:`p`.
   :param int nxi: Number of distribution grid points in the spherical coordinate :math:`\xi`.
   :param float pmax: Value of upper flux grid boundary, expressed in the spherical coordinate :math:`p`.

.. py:method:: set(enabled=True, ttype=MOMENTUMGRID_TYPE_PXI, np=100, nxi=1, pmax=None)

   Set all settings for this momentum grid in one go after creating the object.

   :param bool enabled: Whether or not the grid is to be used during the simulation.
   :param int type: Type of momentum grid (currently, only ``MOMENTUMGRID_TYPE_PXI``, for spherical coordinates, is supported).
   :param int np: Number of distribution grid points in the spherical coordinate :math:`p`.
   :param int nxi: Number of distribution grid points in the spherical coordinate :math:`\xi`.
   :param float pmax: Value of upper flux grid boundary, expressed in the spherical coordinate :math:`p`.

.. py:method:: setEnabled(enabled=True)

   Specifies whether or not this momentum grid should be enabled and used
   during the DREAM simulation.

   :param bool enabled: Whether or not the grid is to be used during the simulation.

.. py:method:: setNp(np)

   Sets the number of points to use for the distribution grid in the spherical
   coordinate :math:`p`.

   :param int np: Number of grid points to use in the spherical coordinate :math:`p`.

.. py:method:: setNxi(nxi)

   Sets the number of points to use for the distribution grid in the spherical
   coordinate :math:`\xi`.

   :param int nxi: Number of grid points to use in the spherical coordinate :math:`\xi`.

.. py:method:: setPmax(pmax)

   Set the value of the upper boundary in the spherical momentum coordinate
   :math:`p`. The value is assigned to the last point on the momentum flux grid.

   :param float pmax: Value of the last momentum flux grid point.

.. py:attribute:: name
.. py:attribute:: pgrid
.. py:attribute:: type
.. py:attribute:: xigrid

Examples
--------
