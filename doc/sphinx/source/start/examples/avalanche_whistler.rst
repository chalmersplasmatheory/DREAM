.. _example-avalanche-whistler:

Quasilinear diffusion from whistler waves
==========================================

This example demonstrates how to set up a DREAM simulation with quasilinear
diffusion due to whistler wave--particle interactions. It is based on the
``examples/avalanche/`` example (``generate_with_fre.py``), but adds
wave-particle interactions using parameters from a QUADRE simulation
(476 MHz whistler wave).

The script configures a kinetic simulation where:

* A uniform electric field accelerates electrons, producing runaway electrons.
* A whistler wave spectrum (476 MHz) is injected during the simulation,
  resonating with high-energy runaway electrons (p ≈ 16, about 7.7 MeV).
* Quasilinear diffusion from the wave scatters the runaway electrons in
  momentum space, potentially suppressing runaway growth.

Physics background
------------------
Whistler waves (also called "chorus waves") are electromagnetic waves that
propagate in magnetised plasmas at frequencies between the ion cyclotron
frequency and the electron cyclotron frequency. In a tokamak disruption, these
waves can be excited by the runaway electron beam itself, and in turn scatter
the runaways through resonant wave--particle interactions.

The resonance condition for whistler waves and electrons is approximately

.. math::

   \omega - k_\parallel v_\parallel = - n \Omega_{\rm ce} / \gamma,

where :math:`\omega` is the wave frequency, :math:`k_\parallel` is the
parallel wavenumber, :math:`v_\parallel` is the parallel electron velocity,
:math:`\Omega_{\rm ce}` is the electron cyclotron frequency, :math:`\gamma` is
the relativistic factor, and :math:`n = -2,-1,0,+1,+2` are the harmonic numbers.

Running the example
-------------------
The script is located at

::

   examples/avalanche_whistler/scripts/generate_with_fre_whistler.py

Run with default parameters:

.. code-block:: bash

   $ python3 examples/avalanche_whistler/scripts/generate_with_fre_whistler.py

Customise parameters:

.. code-block:: bash

   $ python3 examples/avalanche_whistler/scripts/generate_with_fre_whistler.py \\
       --amplitude 1e3 --a 0.3 --R 1.67 --E 0.05

View all available options:

.. code-block:: bash

   $ python3 examples/avalanche_whistler/scripts/generate_with_fre_whistler.py --help

The script produces two files:

* ``quasilinear_whistler_settings.h5`` — the DREAM settings file (can be reused
  with ``dreami`` directly)
* The output HDF5 file (default: ``../outputs/quasilinear_whistler_output.h5``)

Command line parameters
-----------------------
.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``--amplitude``
     - ``1e3``
     - Wave amplitude in normalised units (δB/:math:`B_0` scaling)
   * - ``--E``
     - ``0.05``
     - Electric field strength [V/m]
   * - ``--n``
     - ``5e18``
     - Electron density [m⁻³]
   * - ``--T``
     - ``2165``
     - Electron temperature [eV]
   * - ``--a``
     - ``0.3``
     - Minor radius [m]
   * - ``--R``
     - ``1.67``
     - Major radius [m]
   * - ``--B0``
     - ``1.4``
     - On-axis magnetic field [T]
   * - ``--Np-hot``
     - ``100``
     - Hot-tail momentum grid points
   * - ``--Np-re``
     - ``200``
     - Runaway momentum grid points
   * - ``--Nxi``
     - ``40``
     - Pitch grid points
   * - ``--tMax``
     - ``2.5``
     - Simulation time [s]
   * - ``--Nt``
     - ``2500``
     - Number of time steps
   * - ``--output``
     - ``../outputs/quasilinear_whistler_output.h5``
     - Output file path
   * - ``--start-inject-time``
     - ``1.0``
     - Wave injection start time [s]
   * - ``--inject-cycle-duration``
     - ``0.2``
     - Wave injection ON+OFF cycle duration [s]
   * - ``--ramp-time``
     - ``0.02``
     - Ramp-up time for wave injection [s]
   * - ``--source``
     - ``off``
     - Enable (``on``) or disable (``off``) kinetic avalanche source

Key configuration details
-------------------------
**Plasma and grid setup:**

* Uniform electric field and temperature
* Fully ionised deuterium plasma
* Hot-tail grid with p_max = 1 and runaway grid with p_max = 50
* 1 radial point (0D-like, uniform plasma)
* Trapped-passing boundary layer grid with dξ_max = 0.05

**Runaway generation:**

* Dreicer generation is disabled (the kinetic simulation naturally captures
  the hot-tail runaway mechanism)
* Avalanche generation (kinetic model) can be toggled via ``--source``:
  when enabled, the avalanche source term is included in the kinetic
  equation for :math:`f_{\rm re}`

**Wave parameters (QUADRE):**

The whistler wave spectrum is defined by the following parameters,
obtained from a QUADRE simulation:

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Parameter
     - Value
     - Description
   * - Frequency
     - 476 MHz
     - Whistler wave frequency
   * - k_total
     - 54.58 m⁻¹
     - Total wavenumber
   * - :math:`k_\parallel`
     - -41 m⁻¹
     - Parallel wavenumber (backward propagation)
   * - k range
     - [51.61, 59.55] m⁻¹
     - Wavenumber spectrum range
   * - θ range
     - [2.38, 2.47] rad
     - Propagation angle range

The quasilinear diffusion operator includes five harmonics
(:math:`n = -2, -1, 0, +1, +2`) and uses a uniform grid of :math:`8 \times 20`
modes in (k, θ) space.

**Radial grid:**

Uses an analytic toroidal geometry (:code:`TYPE_ANALYTIC_TOROIDAL`) with

.. code-block:: text

   ψ = 0.0001     (inverse aspect ratio shaping)
   G/R₀ = B₀      (poloidal current function)

**Solver:**

* Linear implicit time integration
* Preconditioner enabled
* Quantities ``nu_s``, ``nu_D``, ``lnLambda`` saved for analysis

Expected output
---------------
The simulation evolves the electron distribution function from the initial
Maxwellian through the acceleration phase and into the runaway regime. When
the whistler wave is injected (at ``start_inject_time``), quasilinear
diffusion should modify the high-energy tail of the distribution and may
reduce the net runaway production rate. The output file contains the full
kinetic distribution functions :math:`f_{\rm hot}` and :math:`f_{\rm re}`,
fluid quantities, and the diagnostic quantities requested via
:code:`ds.other.include(...)`.

Plotting the runaway current over time with and without wave injection can
illustrate the effect of the wave--particle interactions.
