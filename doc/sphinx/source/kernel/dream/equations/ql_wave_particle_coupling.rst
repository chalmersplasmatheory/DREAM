.. _ql-wave-particle-coupling:

WaveParticleCoupling
====================

The :cpp:class:`WaveParticleCoupling` class calculates the coupling strength
:math:`|\Psi_{n,k}|^2` between the wave mode :math:`(k, \theta_k)` and
particles at a given phase-space point :math:`(p, \xi)` for harmonic
:math:`n`.

.. code-block:: cpp

   class WaveParticleCoupling {
   public:
       WaveParticleCoupling(real_t omega_pe, real_t omega_ce, real_t n_e);
       
       real_t calculateCouplingStrength(
           real_t p, real_t xi, real_t k, real_t theta_k, int n,
           const WhistlerDispersion &dispersion
       ) const;
   };

Coupling strength formula
-------------------------
The coupling strength follows the QUADRE implementation (``inject.py``,
based on quasi-linear theory). For a resonant wave-particle interaction:

.. math::

   |\Psi_{n,k}|^2 =
   \frac{1}{|\partial\omega/\partial k - v_\parallel|}
   \cdot \frac{k^2 \omega_{\rm pe}^2}{n_e \pi}
   \cdot W^2,

where the weight :math:`W` encodes the wave polarisation and finite
Larmor radius effects:

.. math::

   W = \frac{n J_n(k_\perp\rho)}{k_\perp\rho}
   + E_z J_n(k_\perp\rho) \frac{\xi}{\sqrt{1-\xi^2}}
   - E_y \frac{J_{n+1}(k_\perp\rho) - J_{n-1}(k_\perp\rho)}{2}.

The terms are:

* :math:`J_n(x)` — Bessel function of the first kind, accounting for
  finite Larmor radius effects
* :math:`\rho = p_\perp / (m_e \Omega_{\rm ce})` — electron gyroradius
* :math:`E_x, E_y, E_z` — wave polarisation vectors from the
  dispersion relation
* :math:`k_\perp = k \sin(\theta_k)` — perpendicular wavenumber
* :math:`\partial\omega/\partial k` — group velocity

Amplitude scaling
-----------------
The coupling strength is scaled by the squared wave amplitude:

.. math::

   |\Psi_{n,k}|^2_{\rm total} = A^2 \cdot |\Psi_{n,k}|^2,

where :math:`A` is the (possibly time-dependent) wave amplitude
in normalised units.

Helper methods
--------------
**Bessel function:**

.. code-block:: cpp

   real_t besselJ(int n, real_t x) const;

Evaluates the Bessel function of the first kind :math:`J_n(x)`, used
in the finite-Larmor-radius coupling terms.

**Group velocity:**

.. code-block:: cpp

   real_t calculateGroupVelocity(real_t k, real_t theta_k,
                                 const WhistlerDispersion &dispersion,
                                 real_t dk = 0.01) const;

Computes :math:`\partial\omega/\partial k` via centred finite differences.

**Gyroradius:**

.. code-block:: cpp

   real_t calculateGyroradius(real_t p, real_t xi) const;

Computes the gyroradius :math:`\rho = p_\perp / (m_e \Omega_{\rm ce})`.

Location
--------
* Header: ``include/DREAM/Equations/Kinetic/WaveParticleCoupling.hpp``
* Source: ``src/Equations/Kinetic/WaveParticleCoupling.cpp``
