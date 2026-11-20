.. _ql-diffusion:

Quasilinear diffusion from wave-particle interactions
======================================================

DREAM implements a quasi-linear diffusion operator for simulating wave-particle
interactions, primarily targeting whistler wave--runaway electron resonance in
tokamak disruptions. The implementation follows the quasi-linear theory
framework and is based on the QUADRE code (Xie 2014).

Physics motivation
------------------
Runaway electrons generated during tokamak disruptions can reach MeV energies
and potentially damage plasma-facing components. Experimental observations have
shown that helicon/whistler wave injection can suppress the growth of runaway
electrons through resonant wave-particle interactions. The key mechanism is:

1. Whistler waves resonate with high-energy runaway electrons, transferring
   parallel momentum to perpendicular momentum (pitch-angle scattering).
2. Higher pitch angle (:math:`\xi`) enhances synchrotron radiation damping,
   causing the electrons to lose energy faster.
3. Electrons scattered into the trapped region (where :math:`|\xi| < \xi_T`)
   are confined by magnetic mirrors and can no longer accelerate freely,
   effectively suppressing runaway generation.

Resonance condition
-------------------
The resonance condition for whistler waves and electrons is:

.. math::

   \omega - k_\parallel v_\parallel = n \frac{\Omega_{\rm ce}}{\gamma},

where

* :math:`\omega` — wave frequency (determined from the dispersion relation)
* :math:`k` — wavenumber magnitude
* :math:`\theta_k` — wave propagation angle with respect to :math:`\mathbf{B}_0`
* :math:`k_\parallel = k \cos(\theta_k)` — parallel wavenumber
* :math:`v_\parallel = (p/\gamma) \xi c` — parallel electron velocity
* :math:`n` — harmonic number (:math:`-2, -1, 0, +1, +2`)
* :math:`p` — normalised momentum (:math:`p = \gamma v / c`)
* :math:`\xi = v_\parallel / v` — pitch-angle cosine
* :math:`\gamma = \sqrt{1 + p^2}` — Lorentz factor
* :math:`\Omega_{\rm ce}` — electron cyclotron frequency

Quasi-linear diffusion operator
--------------------------------
The quasi-linear diffusion term enters the bounce-averaged drift-kinetic
equation as an additional diffusion operator:

.. math::

   \left.\frac{\partial f}{\partial t}\right|_{\rm QL} =
   \frac{1}{p^2} \frac{\partial}{\partial p}
   \left[ p^2 \left( D_{pp} \frac{\partial f}{\partial p}
   + D_{p\xi} \frac{\partial f}{\partial \xi} \right) \right]
   + \frac{\partial}{\partial \xi}
   \left[ \left( D_{p\xi} \frac{\partial f}{\partial p}
   + D_{\xi\xi} \frac{\partial f}{\partial \xi} \right) \right].

The diffusion coefficients are computed as sums over all resonant wave modes
and harmonic numbers:

.. math::

   D_{pp}  &= \sum_{n,k} |\Psi_{n,k}|^2, \\
   D_{p\xi} &= \sum_{n,k} |\Psi_{n,k}|^2 \cdot \left(\frac{k_\parallel v_\perp}{\omega}\right), \\
   D_{\xi\xi} &= \sum_{n,k} |\Psi_{n,k}|^2 \cdot \left(\frac{k_\parallel v_\perp}{\omega}\right)^2,

where :math:`|\Psi_{n,k}|^2` is the wave-particle coupling strength.

Implementation architecture
---------------------------
The quasi-linear diffusion implementation consists of four main classes,
which work together in the :cpp:class:`QuasilinearDiffusionTerm`:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Class
     - Responsibility
     - Source file
   * - :ref:`ql-resonance-solver`
     - Solves the resonance condition to find exact resonant momenta
       for given wave parameters and pitch angle
     - ``ResonanceSolver.cpp``
   * - :ref:`ql-dispersion`
     - Computes wave frequency :math:`\omega(k, \theta_k)` from the
       whistler dispersion relation; provides polarisation vectors
     - ``WhistlerDispersion.cpp``
   * - :ref:`ql-wave-particle-coupling`
     - Calculates the coupling strength :math:`|\Psi_{n,k}|^2` using
       Bessel functions and wave polarisation
     - ``WaveParticleCoupling.cpp``
   * - :ref:`ql-diffusion-coefficients`
     - Assembles the :math:`D_{pp}, D_{p\xi}, D_{\xi\xi}` coefficients
       on the FVM flux grid with REVERSE method
     - ``QuasilinearDiffusionTerm.cpp``

Additionally, the :ref:`ql-builder` module handles construction of the
``QuasilinearDiffusionTerm`` from settings, supporting both on-the-fly
computation and pre-computed matrix (HDF5) modes.

.. toctree::
   :maxdepth: 2
   :hidden:

   ql_resonance_solver
   ql_dispersion
   ql_wave_particle_coupling
   ql_diffusion_coefficients
   ql_builder

Key features
------------
* **REVERSE resonance solving**: instead of iterating over grid points and
  testing for resonance, the method solves for the exact resonant momentum
  for each wave mode, ensuring grid-independent results.
* **Multiple harmonics**: supports :math:`n = -2, -1, 0, +1, +2`.
* **Two modes of operation**: on-the-fly computation (full physics) or
  pre-computed matrix (fast, from HDF5).
* **Periodic wave injection**: time-dependent amplitude with ramp-up,
  duty cycle control and start delay.
* **NaN/Inf validation**: automatic detection of invalid diffusion
  coefficients with detailed diagnostic output.
* **FVM flux-grid interpolation**: triangular weighting distributes
  resonant contributions smoothly across neighbouring grid cells.
