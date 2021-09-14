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

.. todo::

   Describe the need to explicitly disabled undesired momentum grids.

.. todo::

   Describe the use of the bi-uniform momentum grids.

.. todo::

   Describe the use of the custom momentum grids.

.. _ds-momentumgrid-trapped:

Trapped/passing boundary
------------------------
In toroidal geometry, particles with
:math:`|\xi_0| < \sqrt{1-B_{\rm min}/B_{\rm max}}` will be trapped and bounce
back and forth in the magnetic field. Here, :math:`B_{\rm max}` is the maximum
magnetic field strength experienced along the flux surface/orbit, and
:math:`B_{\rm min}` is the minimum magnetic field strength, at which point the
particle instantaneously has :math:`\xi = \xi_0`.

To accurately resolve a distribution function in a toroidal magnetic, the set of
points :math:`\xi_{0,{\rm T}}` referred to as the trapped-passing
boundary---which separates the particles which are trapped from those that are
passing in momentum space---must be resolved very accurately (i.e. we must place
grid points very close to these points). DREAM provides a particular grid type
for automatically locating the trapped-passing boundary and placing grid points
in appropriate locations. To use this grid, call the method
:py:meth:`DREAM.Settings.MomentumGrid.MomentumGrid.setTrappedPassingBoundaryLayerGrid`
on the desired momentum grid.

Since the grid is automatically assembled based on the location of the
trapped-passing boundaries, there are limited options for the user to customize
the grid spacing. The method 
:py:meth:`DREAM.Settings.MomentumGrid.MomentumGrid.setTrappedPassingBoundaryLayerGrid`
takes four arguments, and typically the user may at least want to specify the
parameter ``dxiMax`` to indicate to DREAM the maximum allowed size of each grid
cell. The user can also fine tune the spacing in the fully-trapped and
fully-passing regions separately using the ``nxiPass`` and ``nxiTrap``
parameters, which set (half) the number of grid points to place in each of the
two regions. Finally, the parameter ``boundaryLayerWidth`` indicates how close
points should be placed to the trapped-passing boundary points
:math:`\xi_{0,{\rm T}}` in order to resolve them accurately, and should
generally be a very small number. If the boundary layer width is too large,
the solution can behave weirdly close to the trapped-passing boundary.

.. note::

   The original version of the method
   :py:meth:`DREAM.Settings.MomentumGrid.MomentumGrid.setTrappedPassingBoundaryLayerGrid`
   required the user to manually specify the location of the trapped-passing
   boundary at every radius using the ``xi0Trapped`` parameter. In more recent
   versions, however, the trapped-passing boundary can be automatically
   calculated by DREAM during initialization. This is practically always the
   desired behaviour, in which case ``xi0Trapped`` need not be specified.

Example
*******
.. code-block:: python

   ds = DREAMSettings()
   ...
   ds.hottailgrid.setTrappedPassingBoundaryLayerGrid(dxiMax=1e-3, boundaryLayerWidth=1e-4)

Object documentation
--------------------

.. autoclass:: DREAM.Settings.MomentumGrid.MomentumGrid
   :members:
   :undoc-members:
   :special-members: __init__, __contains__, __getitem__

Examples
--------
Run with a basic hot-tail grid and no fluid runaways:

.. code-block:: python

   from DREAM.DREAMSettings import DREAMSettings
   import DREAM.Settings.CollisionHandler as Collisions

   ds = DREAMSettings()
   ...
   ds.hottailgrid.set(np=500, nxi=10, pmax=5)
   ds.runawaygrid.setEnabled(False)
   ...


