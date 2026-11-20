.. _ql-resonance-solver:

ResonanceSolver
===============

The :cpp:class:`ResonanceSolver` class solves the wave-particle resonance
condition to find exact resonant momenta for given wave parameters and pitch
angle. It implements the **REVERSE** approach: instead of checking whether a
grid point is resonant (which can miss resonances that fall between grid
nodes), it solves for the exact resonant :math:`p` given
:math:`(k, \theta_k, n, \xi)`.

.. code-block:: cpp

   class ResonanceSolver {
   public:
       ResonanceSolver(real_t omega_ce_val);
       
       std::vector<real_t> findResonantP(
           real_t k, real_t theta_k, int n, real_t xi,
           const WhistlerDispersion &dispersion,
           real_t p_min = 0.0, real_t p_max = 100.0
       ) const;
   };

Resonance condition
-------------------
The solver finds momenta :math:`p` satisfying the resonance condition:

.. math::

   \omega(k, \theta_k) - k_\parallel v_\parallel - n \frac{\Omega_{\rm ce}}{\gamma} = 0,

where :math:`v_\parallel = (p/\gamma)\,\xi\,c` and
:math:`k_\parallel = k \cos(\theta_k)`.

To avoid division by :math:`\gamma` (which can cause numerical instability
at low momenta) and to eliminate spurious negative roots from squaring,
the solver uses the equivalent residual function:

.. math::

   f(p) = \omega \gamma(p) - k_\parallel c p \xi - n \Omega_{\rm ce} = 0.

Root-finding algorithm
----------------------
``findResonantP()`` proceeds in three steps:

1. **Sample the residual** — evaluate :math:`f(p)` on a uniform grid of
   200 sample points from :math:`p = 0.01` to :math:`p = 100`.
   Detect sign changes (:math:`f(p_i) \cdot f(p_{i+1}) < 0`) which
   bracket roots.

2. **Bisection** — for each bracketed interval, apply the bisection method
   with up to 50 iterations. Convergence criteria:

   * :math:`|f(p)| < 10^{-6}`
   * Relative change in :math:`p < 10^{-8}`

3. **Verification** — verify the candidate root against the original
   resonance condition.

   :math:`\omega - k_\parallel v_\parallel + n \Omega_{\rm ce} / \gamma = 0`.

The function returns a vector of all resonant momenta found (may be empty
if no resonance exists for the given parameters).

Location
--------
* Header: ``include/DREAM/Equations/Kinetic/ResonanceSolver.hpp``
* Source: ``src/Equations/Kinetic/ResonanceSolver.cpp``
