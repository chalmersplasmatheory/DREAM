.. _ds-eqsys-E_field:

ElectricField
=============
The ``ElectricField`` class holds settings for the parallel electric field
:math:`E_\parallel` solved for by DREAM. The electric field can be solved for
in two different ways:

(1) By prescribing the electric field profile in time :math:`E_\parallel = \tilde{E}(t,r)` (``TYPE_PRESCRIBED``)
(2) By solving the induction equation (``TYPE_SELFCONSISTENT``):

   .. math::

       \begin{cases}
           \frac{\partial\psi_{\rm p}}{\partial t} &= V_{\rm loop},\\
           \mu_0\frac{j_\parallel}{B}\left\langle \boldsymbol{B}\cdot\nabla\varphi \right\rangle &=
           \frac{1}{V'}\frac{\partial}{\partial r}\left[
               V'\left\langle \frac{\left|\nabla r\right|^2}{R^2} \right\rangle
               \frac{\partial\psi_{\rm p}}{\partial r}
           \right]
       \end{cases}

   where :math:`\psi_{\rm p}` is the poloidal flux, :math:`V_{\rm loop}\equiv 2\pi RE_\parallel`
   the loop voltage (at major radius :math:`R`), :math:`\mu_0` is the vacuum permeability,
   :math:`\boldsymbol{B}` the magnetic field, :math:`\varphi` the toroidal angle, :math:`V'` the
   spatial Jacobian, and :math:`r` the minor radius in the outer midplane. Angle brackets
   :math:`\langle\cdot\rangle` denote a flux-surface average.
