.. _ql-diffusion-coefficients:

Diffusion coefficient calculation (REVERSE method)
==================================================

The core of the quasi-linear diffusion implementation is the
``calculateDiffusionCoefficientsReverse()`` method in
:cpp:class:`QuasilinearDiffusionTerm`, which assembles the diffusion
coefficients :math:`D_{pp}, D_{p\xi}, D_{\xi\xi}` on the FVM flux grid.

Algorithm overview
------------------
Instead of iterating over grid points and testing whether each point is
resonant (a "forward" approach that can miss resonances between grid
nodes), the **REVERSE** method works as follows:

1. For each wave mode :math:`(k, \theta_k)`, iterate over pitch-angle
   grid points :math:`\xi_j`.
2. For each harmonic number :math:`n \in \{-2,-1,0,+1,+2\}`, call
   :cpp:func:`ResonanceSolver::findResonantP` to obtain the exact
   resonant momentum :math:`p_{\rm res}`.
3. Locate the FVM grid cell containing :math:`p_{\rm res}` by binary
   search on the flux-face grid.
4. Distribute the resonant contribution to neighbouring grid nodes
   using triangular weights.
5. Accumulate to :math:`D_{pp}, D_{p\xi}, D_{\xi\xi}` on the
   appropriate FVM flux grids (f1 for :math:`D_{pp}`, f2 for
   :math:`D_{p\xi}` and :math:`D_{\xi\xi}`).

Triangular weighting
--------------------
To smoothly distribute the resonant contribution without introducing
grid-dependent artefacts, the contribution of each resonance is spread
over the neighbouring grid cells with a triangular weight function:

.. math::

   w &= 0.5 \, \Delta p, \\
   w_{\rm low}  &= \max\left(0, 1 - \frac{d_{\rm low}}{w}\right), \\
   w_{\rm high} &= \max\left(0, 1 - \frac{d_{\rm high}}{w}\right),

where :math:`d_{\rm low} = |p_{\rm res} - p_{\rm low}|` and
:math:`d_{\rm high} = |p_{\rm res} - p_{\rm high}|` are distances
to the bounding flux-face nodes. The weights are normalised so that
:math:`w_{\rm low} + w_{\rm high} = 1`.

Geometric factors
-----------------
At each resonant point, the coupling strength :math:`|\Psi_{n,k}|^2` is
multiplied by geometric factors to obtain the three diffusion coefficients:

.. math::

   D_{pp}       &\propto |\Psi_{n,k}|^2 \cdot (1 - \xi^2), \\
   D_{p\xi}     &\propto -|\Psi_{n,k}|^2 \cdot \sqrt{1-\xi^2} \cdot
                  \left(\xi - \frac{k_\parallel v_\parallel}{\omega}\right), \\
   D_{\xi\xi}   &\propto  |\Psi_{n,k}|^2 \cdot
                  \left(\xi - \frac{k_\parallel v_\parallel}{\omega}\right)^2.

The factor :math:`(1-\xi^2)` for :math:`D_{pp}` follows the form derived
in Zehua Guo (2024) Eq. 11.

Quadrature normalisation
-------------------------
The wave spectrum integral :math:`\int dk \, d(\cos\theta_k)` is
implemented as a discrete sum over a uniform grid of :math:`N_k \times
N_{\theta}` modes:

.. math::

   \iint dk \, d(\cos\theta_k) \approx \sum_{m} \Delta k \cdot \Delta(\cos\theta_k),

where :math:`\Delta k` and :math:`\Delta(\cos\theta_k)` are the grid
spacings in wavenumber and propagation angle, respectively. Each mode's
contribution is scaled by its cell area in :math:`(k, \cos\theta_k)` space.

Caching
-------
The diffusion coefficients depend only on the wave spectrum and plasma
parameters, **not** on the distribution function :math:`f(p, \xi)`. They
are therefore computed once and cached. The cache is invalidated when:

* The wave amplitude changes (e.g. due to periodic injection timing),
* The grid dimensions change (on mesh refinement),
* The external amplitude is modified through the settings interface.

Periodic wave injection
-----------------------
The time-dependent amplitude follows a periodic schedule:

.. code-block:: cpp

   real_t calculateEffectiveAmplitude(real_t t) const;

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Behaviour
   * - ``start_inject_time``
     - Injection starts at this time (if :math:`\ge 0`); before that, amplitude = 0
   * - ``ramp_time``
     - Linear ramp-up from 0 to full amplitude over this duration at each ON cycle
   * - ``inject_cycle_duration``
     - Total ON+OFF cycle length; 50% duty cycle (ON = first half, OFF = second half)
   * - ``inject_cycle_duration <= 0``
     - Continuous injection (no periodicity) after start

FVM matrix assembly
-------------------
The cached coefficients are interpolated onto the FVM diffusion matrices
each time step:

* :math:`D_{pp}` is stored on the **f1** (p-flux) grid:
  :math:`(N_p+1) \times N_\xi`
* :math:`D_{p\xi}` and :math:`D_{\xi\xi}` are stored on the **f2**
  (:math:`\xi`-flux) grid:
  :math:`N_p \times (N_\xi+1)`

These matrices enter the FVM discretisation as the :math:`D_{11}`,
:math:`D_{12}` and :math:`D_{22}` components of the
:cpp:class:`FVM::DiffusionTerm`.

NaN/Inf validation
------------------
After assembly, the coefficients are automatically scanned for NaN and Inf
values. If invalid entries are detected, a detailed diagnostic message is
printed, including the grid indices and the suspected cause.

Location
--------
* Header: ``include/DREAM/Equations/Kinetic/QuasilinearDiffusionTerm.hpp``
* Source: ``src/Equations/Kinetic/QuasilinearDiffusionTerm.cpp``
