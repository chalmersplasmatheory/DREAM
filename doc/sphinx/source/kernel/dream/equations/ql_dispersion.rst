.. _ql-dispersion:

WhistlerDispersion
==================

The :cpp:class:`WhistlerDispersion` class computes the wave frequency
:math:`\omega(k, \theta_k)` from the whistler wave dispersion relation.
It also provides the polarisation vectors and normalisation factors
needed by :cpp:class:`WaveParticleCoupling`.

.. code-block:: cpp

   class WhistlerDispersion {
   public:
       WhistlerDispersion(real_t B0, real_t density,
                         real_t ion_mass_ratio = 1836.0,
                         real_t Zeff = 1.0,
                         bool use_simple = false);
       
       real_t calculateOmega(real_t k, real_t theta) const;
       real_t calculateOmegaSimple(real_t k, real_t theta) const;
       real_t calculateKperpSimple(real_t omega, real_t k_par) const;
       real_t calculateGroupVelocity(real_t k, real_t theta) const;
       
       void calculatePolarization(real_t omega, real_t k, real_t theta,
                                  real_t &Ex, real_t &Ey, real_t &Ez) const;
       real_t calculateDenominator(real_t omega, real_t k, real_t theta) const;
       
       // Accessors
       real_t getOmegaCE() const;
       real_t getOmegaPE() const;
       real_t getB0() const;
   };

Simplified whistler dispersion
-------------------------------
For frequencies satisfying :math:`\Omega_{\rm ci} \ll \omega \ll \Omega_{\rm ce}`,
the cold plasma whistler dispersion simplifies to:

.. math::

   \omega = k |k_\parallel| w, \qquad
   w = \frac{\Omega_{\rm ce} c^2}{\omega_{\rm pe}^2},

where :math:`\omega_{\rm pe}` is the electron plasma frequency. This
approximation is fast and avoids the full eigenmode solve; it is enabled
by setting ``use_simple = true`` in the constructor.

Cold plasma dispersion (Stix formalism)
----------------------------------------
When the simplified relation is not used, the dispersion is solved via the
cold plasma dispersion relation in Stix notation:

.. math::

   A n^4 - B n^2 + C = 0,

where :math:`n = kc/\omega` is the refractive index and the coefficients
:math:`A, B, C` are functions of the Stix parameters:

.. math::

   R &= 1 - \sum_s \frac{\omega_{ps}^2}{\omega(\omega + \Omega_{cs})}, \\
   L &= 1 - \sum_s \frac{\omega_{ps}^2}{\omega(\omega - \Omega_{cs})}, \\
   P &= 1 - \sum_s \frac{\omega_{ps}^2}{\omega^2}, \\
   S &= \frac{R + L}{2}, \quad D = \frac{R - L}{2}.

The sum is over all species :math:`s` (electrons and ions).

PDRF method
-----------
The class also supports solving the full multi-fluid dispersion relation
using the generalised eigenvalue method from Xie's PDRF code (2014). In
this approach, a :math:`(4S+6) \times (4S+6)` matrix is constructed from
the linearised fluid equations and Maxwell's equations:

.. math::

   \mathbf{M} \cdot \mathbf{x} = \omega \cdot \mathbf{A} \cdot \mathbf{x},

where :math:`\mathbf{x}` contains :math:`4S` perturbation variables
(:math:`\delta n, \delta v_x, \delta v_y, \delta v_z`) for each of
:math:`S` species, plus 6 electromagnetic field components
(:math:`E_x, E_y, E_z, B_x, B_y, B_z`). The whistler branch is selected
from the resulting eigenvalues.

Group velocity
--------------
The group velocity is computed via numerical differentiation:

.. math::

   \frac{\partial \omega}{\partial k} \approx
   \frac{\omega(k + \Delta k, \theta_k) - \omega(k - \Delta k, \theta_k)}{2 \Delta k}.

Polarisation
------------
The polarisation vectors :math:`(E_x, E_y, E_z)` are obtained from the
cold plasma dielectric tensor eigenvectors, following the QUADRE convention.

Location
--------
* Header: ``include/DREAM/Equations/Kinetic/WhistlerDispersion.hpp``
* Source: ``src/Equations/Kinetic/WhistlerDispersion.cpp``

References
----------
* Xie, H. S., *PDRF: A general dispersion relation solver for magnetized
  multi-fluid plasma*, Computer Physics Communications **185**, 670-675 (2014).
