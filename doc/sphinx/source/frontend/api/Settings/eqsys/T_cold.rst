.. _ds-eqsys-T_cold:

ColdElectronTemperature
=======================
The :py:class:`DREAM.Settings.Equations.ColdElectronTemperature.ColdElectronTemperature`
holds settings for the (cold) electron temperature. The unknown quantity is
called ``T_cold`` and can be evolved either according to a prescribed function,
or according to an energy balance equation.

The electron temperature and electron thermal energy are related through

.. math::

   W_{\rm cold} = \frac{3}{2}n_{\rm cold}T_{\rm cold}.


.. contents:: Page overview
   :local:
   :depth: 3

Prescribed evolution
--------------------
The temperature can be prescribed in both space and time. The full syntax for
prescribing the electron temperature is

.. code-block::

   ds.eqsys.T_cold.setPrescribedData(temperature=T_cold, radius=T_cold_r, times=T_cold_t)

where ``T_cold`` is an array of shape ``(nt, nr)`` specifying the cold electron
temperature in space and time, ``T_cold_r`` is a vector with ``nr`` elements
representing the radial grid on which the temperature is prescribed, and
``T_cold_t`` is a vector with ``nt`` elements representing the time grid on
which the temperature is prescribed.

For convenience, it is possible to prescribe a constant and radially uniform
temperature by only specifying it as a scalar:

.. code-block::

   ds.eqsys.T_cold.setPrescribedData(1000)

This will cause the temperature to be :math:`T_{\rm cold} = 1\,\mathrm{keV}` at
every radius, in every time step, during the simulation.

Self-consistent evolution
-------------------------
The temperature can be evolved self-consistently, i.e. by calculating the energy
loss of the plasma through various processes. In its most general form, the
energy balance equation solved by DREAM is

.. math::
   :label: eq_energy

   \frac{\partial W_{\rm cold}}{\partial t} &=
       \frac{j_\Omega}{B}\left\langle \boldsymbol{E}\cdot\boldsymbol{B} \right\rangle -
       \left\langle n_{\rm cold} \right\rangle\sum_i\sum_{j=0}^{Z_i-1} n_i^{(j)}L_i^{(j)} +
       \left\langle Q_c \right\rangle +\\
       &+ \frac{1}{V'}\frac{\partial}{\partial r}\left[
           V'\frac{3\left\langle n_{\rm cold} \right\rangle}{2}\left(
               -A_WT_{\rm cold} + D_W\frac{\partial T_{\rm cold}}{\partial r}
           \right)
       \right]

Here, the terms included in the most complete form of the energy balance
equation are

- :math:`\frac{j_\Omega}{B}\langle\boldsymbol{E}\cdot\boldsymbol{B}\rangle`:
  Ohmic heating (energy transferred to the plasma from the electric field).
- :math:`\langle n_{\rm cold} \rangle\sum_i\sum_{j=0}^{Z_i-1} n_i^{(j)} L_i^{(j)}`:
  Energy loss by inelastic atomic processes. The rate coefficient :math:`L_i^{(j)}`
  is given by :math:`L_i^{(j)} = L_{\rm line} + L_{\rm free} +
  \Delta W_i^{(j)}(I_i^{(j)} - R_i^{(j)})`, where :math:`L_{\rm line}` is the
  radiated power by line emission, :math:`L_{\rm free}` is the radiated power
  by recombination emission and bremsstrahlung, :math:`\Delta W_i^{(j)}` is the
  ionization threshold (taken from NIST), and :math:`I_i^{(j)}` and
  :math:`R_i^{(j)}` are the rates of ionization and recombination for the given
  ion charge state. Unless otherwise specified, coefficients are taken from
  OpenADAS.
- :math:`\langle Q_c \rangle`: Collisional heat transfer to the cold electrons,
  from hot and runaway electrons, as well as ions:

  .. math::

     \left\langle Q_c \right\rangle &= \int\mathrm{d}p \int_{-1}^1\mathrm{d}\xi_0
     \frac{\mathcal{V}'}{V'}\Delta\dot{E}_{ee} f_{\rm hot/re} + \sum_i Q_{ei},\\
     Q_{kl} &= \frac{\left\langle nZ^2 \right\rangle_k\left\langle nZ^2 \right\rangle_l
        e^4\ln\Lambda_{kl}}{\left(2\pi\right)^{3/2}\epsilon_0^2 m_k m_l}
        \frac{T_k - T_l}{\left(\frac{T_k}{m_k} + \frac{T_l}{m_l} \right)^{3/2}},\\
     \Delta\dot{E}_ee &= 4\pi n_{\rm cold} r_0^2 \ln\Lambda_{ee}\frac{m_ec^4}{v},
- **Radial transport** of heat. DREAM uses and advection-diffusion heat
  transport model with either arbitrary coefficients :math:`A_W` and
  :math:`D_W`, or a Rechester-Rosenbluth diffusion model requiring the magnetic
  perurbation level :math:`\delta B/B` to be specified.

To solve for the electron temperature self-consistently, the user should call
the method ``setType()`` with the argument ``TYPE_SELFCONSISTENT``.
Additionally, the temperature profile at :math:`t=0` must be specified:

.. code-block:: python

   import DREAM.Settings.Equations.ColdElectronTemperature as Tcold

   ds = DREAMSettings()
   ...
   ds.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
   ds.eqsys.T_cold.setInitialProfile(temperature=T0, radius=T0_r)

If a uniform initial temperature profile is desired, the parameter
``temperature`` in the call to ``setInitialProfile()`` can be a scalar, in which
case the parameter ``radius`` can also be omitted.


Recombination radiation
***********************
By default, radiation due to recombination is neglected in the energy balance
equation :eq:`eq_energy` (i.e. the coefficient :math:`L_{\rm free}` contains
only the contribution from bremsstrahlung, and the coefficient :math:`R_i^{(j)}`
is dropped altogether). To account for radiation due to recombination, the user
must call the method ``setRecombinationRadiation`` with the option
``RECOMBINATION_RADIATION_INCLUDED`` or ``True``.

.. code-block:: python

   import DREAM.Settings.Equations.ColdElectronTemperature as Tcold

   ds = DREAMSettings()
   ...
   ds.eqsys.T_cold.setRecombinationRadiation(Tcold.RECOMBINATION_RADIATION_INCLUDED)
   # or...
   ds.eqsys.T_cold.setRecombinationRadiation(True)


Class documentation
-------------------

.. autoclass:: DREAM.Settings.Equations.ColdElectronTemperature.ColdElectronTemperature
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__
