.. _ds-eqsys-n_re:

RunawayElectrons
================
The ``RunawayElectons`` object, located under ``.eqsys.n_re`` in the
:ref:`DREAMSettings` object is used to specify settings for the fluid quantity
``n_re``, representing the density of runaway electrons in simulations. The
quantity is always present in DREAM simulations and is evolved using the
equation

.. math::

   \frac{\mathrm{d} n_{\rm RE}}{\mathrm{d} t} = \Phi^{(p)}_{\rm hot} +
   \Gamma_{\rm Ava}n_{\rm RE} + \gamma_{\rm Dreicer} + \gamma_{\rm hottail} + 
   \gamma_{\rm tritium} + \gamma_{\rm Compton} +
   \frac{1}{V'}\frac{\partial}{\partial r}\left[
       V'\left( An_{\rm RE} +  D\frac{\partial n_{\rm RE}}{\partial r} \right)
   \right]

where

- :math:`\Phi^{(p)}_{\rm hot}` is the flux of particles from the (kinetic) hot-tail grid.
- :math:`\Gamma_{\rm Ava}` is the avalanche growth rate.
- :math:`\gamma_{\rm Dreicer}` is the runaway rate due to the Dreicer mechanism.
- :math:`\gamma_{\rm hottail}` is the runaway rate due to the hottail mechanism.
- :math:`\gamma_{\rm tritium}` is the runaway rate due to tritium decay.
- :math:`\gamma_{\rm Compton}` is the runaway rate due to Compton scattering.

and the last term represents advective-diffusive radial transport. Each term can
be individually enabled/disabled using the methods described further down on
this page.

.. note::

   The Dreicer generation rate :math:`\gamma_{\rm Dreicer}` should typically
   be disabled in simulations including the hot-tail grid where the thermal population 
   is resolved kinetically. This is because the
   Dreicer generation will then be modelled using a kinetic equation and the
   Dreicer generation rate become a part of the hot electron flux
   :math:`\Phi^{(p)}_{\rm hot}`.


.. note::

   Note that the runaway density :math:`n_{\rm RE}` is **never** defined in
   terms of the runaway electron distribution function :math:`f_{\rm RE}` in
   DREAM. The quantities are instead evolved separately and set up in such a
   way that their local densities should agree in the end if the simulation was
   set up in a consistent way.

   *It is possible---and perfectly okay---for* :math:`n_{\rm RE}` *to differ
   from the density moment of* :math:`f_{\rm RE}` *if the latter loses a
   significant number of particles at* :math:`p=p_{\rm max}`.


.. contents:: Page overview
   :local:
   :depth: 3



Initialization
--------------
By default, the runaway density is set to zero at the beginning of a DREAM
simulation. It is however possible to start the simulation with a non-zero
runaway density profile, and this is done using the ``setInitialProfile()``
method.

.. warning::

   If the runaway grid is enabled, one should also make sure to appropriately
   initialize the runaway electron distribution function when setting a non-zero
   initial runaway density. The density moment of the runaway electron
   distribution function :math:`f_{\rm RE}` should always agree with the runaway
   electron density :math:`n_{\rm RE}`.

Example
^^^^^^^
The following example illustrates how to initialize the runaway electron
density:

.. code-block:: python

   import numpy as np

   ds = DREAMSettings()
   ...
   # Define radial grid for the initial profile
   r = np.linspace(r0, r1, nr)

   # Generate the runaway electron density profile
   nre = np.array([...])        # shape (nr,)

   ds.eqsys.n_re.setInitialProfile(density=nre, r=r)

If a uniform initial density profile is desired one can specify the ``density``
parameter as a scalar value instead:

.. code-block:: python

   import numpy as np

   ds = DREAMSettings()
   ...
   ds.eqsys.n_re.setInitialProfile(density=2e16)


Runaway generation rates
------------------------
DREAM includes a number of common runaway electron source terms for use in both
pure fluid as well as combined fluid-kinetic simulations.

Avalanche
^^^^^^^^^
Runaway electrons can produce new runaway electrons by colliding with thermal
electrons and transferring a sufficent amount of energy to these for them to
become runaway electrons, while the original runaway electron remains in the
runaway region. This runaway mechanism is commonly known as the *avalanche
mechanism* and will lead to an exponential growth in the number of runaway
electrons.

DREAM contains three different models for avalanche generation. The first two
are fluid models, while the third is the kinetic source term originally derived
by `Rosenbluth and Putvinski (1997) <https://doi.org/10.1088/0029-5515/37/10/I03>`_.

.. note::

   Note that also fluid growth rates can be used in simulations involving the
   runaway grid. In these simulations, the runaways produced with the fluid
   source term are added to the :math:`\xi=1` cell on the runaway grid (if
   :math:`E > 0`, otherwise to the :math:`\xi=-1` cell).

   Similarly, the kinetic source term can be used in fluid simulations and is
   then integrated to yield the number of fluid runaways produced.

+----------------------------------+---------------------------------------------------------------------------------------+
| Option                           | Description                                                                           |
+==================================+=======================================================================================+
| ``AVALANCHE_MODE_NEGLECT``       | Do **not** include avalanche generation in the simulation.                            |
+----------------------------------+---------------------------------------------------------------------------------------+
| ``AVALANCHE_MODE_FLUID``         | Modified ``HESSLOW`` model                                                            |
+----------------------------------+---------------------------------------------------------------------------------------+
| ``AVALANCHE_MODE_FLUID_HESSLOW`` | Implements Eq (14) in `Hesslow NF (2019) <https://doi.org/10.1088/1741-4326/ab26c2>`_ |
+----------------------------------+---------------------------------------------------------------------------------------+
| ``AVALANCHE_MODE_KINETIC``       | Kinetic source term derived by Rosenbluth and Putvinski (1997).                       |
+----------------------------------+---------------------------------------------------------------------------------------+

The ``FLUID`` model modifies the ``FLUID_HESSLOW`` model by replacing the denominator 
:math:`\sqrt{4+\bar\nu_s\bar\nu_D}` appearing in (14) of the Hesslow paper by 
:math:`\sqrt{4\bar\nu_s+\bar\nu_s\bar\nu_D}`, which increases the accuracy for 
nearly neutral plasmas dominated by hydrogen collisions.

.. note::
   Running the ``dream_avalanche`` physics test with the --plot flag 
   compares the ``FLUID`` and ``FLUID_HESSLOW`` models with the 
   kinetic calculation. 


Example
*******
The following example illustrates how to enable the kinetic runaway source term:

.. code-block:: python

   import DREAM.Settings.Equations.RunawayElectrons as Runaways

   ds = DREAMSettings()
   ...
   ds.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_KINETIC)


Dreicer
^^^^^^^
The first runaway mechanism to be discovered is the so-called *Dreicer runaway
mechanism* which occurs in a plasma as soon as a sufficiently strong electric
is applied. The electric field immediately accelerates all electrons with
momentum :math:`p > p_{\rm c}` to even higher momentum, leaving a "gap" in the
electron distribution function. Since collisions seek to maintain the Maxwellian
form for the distribution, the electrons soon reorganize to replace the
accelerated electrons. However, the new electrons which diffused to fill in the
gap now find themselves with a momentum :math:`p > p_{\rm c}` and are thus also
immediately accelerated to higher energy. This leads to gradual diffusive leak
of particles into the runaway region.

Since the Dreicer mechanism is naturally obtained as the balance between
collisions and electric field acceleration, the mechanism is automatically and
always present in all kinetic simulations. To enable Dreicer generation in pure
fluid simulations, the options described in the table below are available.

+---------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Option                                | Description                                                                                                                                                                              |
+=======================================+==========================================================================================================================================================================================+
| ``DREICER_RATE_DISABLED``             | Do not include Dreicer generation in the simulation.                                                                                                                                     |
+---------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``DREICER_RATE_CONNOR_HASTIE_NOCORR`` | Use the formula derived by `Connor and Hastie (1975) <https://doi.org/10.1088/0029-5515/15/3/007>`_, excluding the relativistic corrections (valid in the limit :math:`E\gg E_{\rm c}`). |
+---------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``DREICER_RATE_CONNOR_HASTIE``        | Use the full formula derived by `Connor and Hastie (1975) <https://doi.org/10.1088/0029-5515/15/3/007>`_.                                                                                |
+---------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``DREICER_RATE_NEURAL_NETWORK``       | Use the neural network constructed by `Hesslow et al (2019) <https://doi.org/10.1017/S0022377819000874>`_.                                                                               |
+---------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Example
*******
The following example illustrates how to enable the kinetic runaway source term:

.. code-block:: python

   import DREAM.Settings.Equations.RunawayElectrons as Runaways

   ds = DREAMSettings()
   ...
   ds.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)


Tritium
^^^^^^^
Runaway electrons can be generated when tritium decays into helium-3 via the
beta decay process

.. math::

   \mathrm{T} \to \ ^3_2\mathrm{He} + \mathrm{e}^- + \bar{\nu}_{\rm e}.

The corresponding runaway rate is given by

.. math::

   \left( \frac{\mathrm{d} n_{\rm RE}}{\mathrm{d} t} \right)_{\rm T} \approx
   \ln 2 \frac{n_{\rm T}}{\tau_{\rm T}} F_\beta\left( \gamma_{\rm c} \right),

where :math:`n_{\rm T}` is the tritium density, :math:`\tau_{\rm T} = 4800\pm 8`
days is the tritium half-life, and :math:`F_\beta(\gamma_{\rm c})` denotes the
fraction of beta electrons generated with an energy above the critical energy
:math:`\gamma_{\rm c}` for runaway to occur.

Tritium runaway generation is enabled with the ``setTritium()`` method of the
``RunawayElectrons`` object, but it is necessary to also provide a tritium ion
species in the :ref:`ions interface <ds-eqsys-ions-tritium>`. The user must make
sure to specify ``tritium=True`` when adding the tritium ion to the simulation
using the ``addIon()`` method.

.. note::

   It is possible to include multiple tritium populations in the simulation,
   simply by providing ``tritium=True`` when adding each of them.

Example
*******
The following example illustrates how to enable the tritium decay runaway
mechanism in a DREAM simulation:

.. code-block:: python

   import DREAM.Settings.Equations.IonSpecies as Ions

   ds = DREAMSettings()
   ...
   # Include source term in equation for n_re
   ds.eqsys.n_re.setTritium(True)

   # Add tritium ion species to list of ions
   ds.eqsys.ions.addIon('T', Z=1, iontype=Ions.IONS_DYNAMIC, n=2e19, tritium=True)


Compton
^^^^^^^
Runaway generation by Compton scattering is modelled via the Klein-Nishina differential
cross section integrated over the runaway region :math:`p>p_c`, analogous to the tritium 
and avalanche calculations. This integrated cross section is given by Equation (29) in
`Martin-Solis NF (2017) <https://doi.org/10.1088/1741-4326/aa6939>`_.
This total cross section :math:`\sigma(p_c,\,E_\gamma)`
depends on plasma parameters through :math:`p_c`and on the incident-photon energy 
:math:`E_\gamma`. The net runaway rate is consequently obtained as the integral 
of the total cross section over the photon energy spectrum :math:`\Phi(E_\gamma)`
(where we assume that bound and free electrons are equally susceptible to 
Compton scattering):

.. math::
   \left(\frac{\partial n_\mathrm{RE}}{\partial t}\right)_\mathrm{Compton} 
   = n_\mathrm{tot} \int \Phi \sigma \, \mathrm{d}E_\gamma

Following Martin-Solis Eq (24), we model the photon energy spectrum with 

.. math::
   \Phi = \Phi_0 \frac{\exp[ - \exp(-z) - z + 1 ]}{5.8844 m_e c^2} 

where :math:`z = [1.2 + \ln(E_\gamma[MeV])]/0.8`, and 
:math:`\Phi_0 = \int \Phi \,\mathrm{d}E_\gamma` is the total photon gamma flux.
For ITER, according to Martin-Solis, this value is :math:`\Phi_0 \approx 10^{18}\,\mathrm{m}^{-2}\mathrm{s}^{-1}`
during the nuclear phase.

The following settings are used to control the Compton source mode:

+--------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| Option                   | Description                                                                                                                       |
+==========================+===================================================================================================================================+
| ``COMPTON_MODE_NEGLECT`` | Do not include Compton scattering in the simulation.                                                                              |
+--------------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ``COMPTON_MODE_FLUID``   | Use the model described in `Martin-Solis NF (2017) <https://doi.org/10.1088/1741-4326/aa6939>`_, with tuneable total photon flux. |
+--------------------------+-----------------------------------------------------------------------------------------------------------------------------------+

Example
*******

The Compton source can be activated following this example: 

.. code-block:: python 

   import DREAM.Settings.Equations.RunawayElectons as Runaways

   ds = DREAMSettings()

   Phi0 = 1e18 # total photon flux in units of m^-2 s^-1, typical of ITER
   ds.eqsys.n_re.setCompton(compton=Runaways.COMPTON_MODE_FLUID, photonFlux=Phi0)



Hottail
^^^^^^^
A fluid model for hottail generation is implemented based on the Master's thesis 
of `Ida Svenningsson (2020) <https://odr.chalmers.se/handle/20.500.12380/300899>`_.
DREAM mainly utilizes the theory of Section 4.2 in the thesis, using the same 
approximations as used in 'ISOTROPIC' (Nxi=1) kinetic mode. 

The hottail runaway generation rate is given by

.. math::

   \frac{\partial n_{\rm re}}{\partial t} = - 4 \pi p_c^2 \dot{p}_c f_0(p_c)

where :math:`p_c` is the critical runaway momentum, :math:`\dot{p}_c` its time 
rate of change and :math:`f_0` is the angle-averaged electron distribution 
evaluated at the critical momentum. In this hottail model, the critical momentum 
is found from the angle-averaged kinetic equation in the Lorentz limit of strong 
pitch-angle scattering:

.. math::

   0 = \frac{1}{3}\left(\frac{E}{E_c}\right)^2 \frac{1}{1+Z_{\rm eff}} \frac{p^3}{\gamma} \frac{\partial f_0}{\partial p} + \frac{\gamma^2}{p^2}f_0,

which for a given distribution :math:`f_0` is solved numerically for the critical 
momentum :math:`p=p_c`. Here we assume that only pitch-angle scattering, electric
field acceleration and collisional slowing down contributes to the electron dynamics.
We have taken the non-relativistic limit of the equation, and neglected screening
corrections which requires a high effective charge in order to be accurate (typically
satisfied during the thermal quench with high-Z material injection). 

The distribution function :math:`f_0` is taken as the analytical hot electron 
distribution described in :ref:`HotElectronDistribution<ds-eqsys-fhot>`.

The following hot tail settings are supported:

+----------------------------------+------------------------------------------------------------------+
| Option                           | Description                                                      |
+==================================+==================================================================+
| ``HOTTAIL_MODE_DISABLED``        | Do **not** include hottail generation in the simulation.         |
+----------------------------------+------------------------------------------------------------------+
| ``HOTTAIL_MODE_ANALYTIC``        | TODO (the low-Z limit, Section 4.1.3 in Svenningsson's thesis)   |
+----------------------------------+------------------------------------------------------------------+
| ``HOTTAIL_MODE_ANALYTIC_ALT_PC`` | The high-Z limit from Section 4.2 in Svenningsson's thesis       |
+----------------------------------+------------------------------------------------------------------+


Example
*******
The hottail generation can be activated if and only if ``f_hot`` is in ``analytical`` mode.

.. code-block:: python 

   import DREAM.Settings.Equations.RunawayElectons as Runaways

   ds = DREAMSettings()

   # Set f_hot mode to analytical and provide initial profiles 
   ds.hottailgrid.setEnabled(False)
   # rn, n0, rT, T0 = ...  get profiles of density and temperature
   ds.eqsys.f_hot.setInitialProfiles(rn0=rn, n0=n0, rT0=rT, T0=T0)

   ds.eqsys.n_re.setHottail(Runaways.HOTTAIL_MODE_ANALYTIC_ALT_PC)



Effective critical electric field
-----------------------------------
DREAM allows the effective critical electric field to be calculated following
the method outlined in `Hesslow PPCF (2018) <https://doi.org/10.1088/1361-6587/aac33e>`_,
which has here been generalized to the case of tokamak geometry 
(i.e. inhomogeneous magnetic fields). The algorithm DREAM implements is relatively
involved, and therefore we briefly describe the method here.


Theory
^^^^^^
Following `Hesslow PPCF (2018) <https://doi.org/10.1088/1361-6587/aac33e>`_ an analytical
expression for the runaway pitch distribution can be obtained under the assumption that 
energy fluxes are negligible compared with pitch fluxes, which is expected to be valid 
near the maximum runaway momentum where energy losses balance the electric field acceleration.

This gives

.. math::

   f_\mathrm{RE} \propto \exp \left[ - \int_{\xi_0}^1 \frac{\{A_{\xi_0}\}}{\{ D_{\xi_0\xi_0} \} } \mathrm{d}\xi_0' \right]

where it is assumed that the dominant bounce averaged terms are the pitch component of the electric field term :math:`A_{\xi_0}` and 
the pitch-angle scattering coefficient :math:`D_{\xi_0\xi_0}`. The exponent can be written in the familar form 

.. math::

   \frac{\{A_{\xi_0}\}}{\{ D_{\xi_0\xi_0} \} } = A g(\xi_0)

where in the cylindrical theory :math:`A = 2eE/(p\nu_D)` and :math:`g = 1-\xi_0`, but in the general case we have

.. math::
   
   A = \frac{2 e \left\langle \boldsymbol{E} \cdot \boldsymbol{B} \right\rangle }{ p \nu_D B_\mathrm{min} }

and :math:`1-\xi_0` is replaced by the function

.. math::

   g = \begin{cases}
      H(\xi_0,\,1), & \xi_T < \xi_0 \leq 1 \\
      H(\xi_T,\,1), & -\xi_T \leq \xi_0 \leq \xi_T \\
      H(\xi_T,\,1) + H(\xi_0, -\xi_T), & -1 \leq \xi_0 < -\xi_T
   \end{cases}

where we introduce the auxiliary function

.. math::

   H(\xi_1,\,\xi_2) = \int_{\xi_1}^{\xi_2} \frac{\xi_0}{\langle \xi \rangle} \mathrm{d}\xi_0

By integrating the kinetic equation over :math:`\mathrm{d}\xi_0 \mathcal{V}'/V'`, it takes the form 

.. math::

   \frac{\partial F_0}{\partial t} + \frac{\partial (U(p)F_0)}{\partial p} = 0

where the distribution-bounce averaged momentum flow (i.e. net acceleration) :math:`U(p)` is 

.. math::

   U = \frac{ \int \mathcal{V}' \{A_p\} f_\mathrm{RE} \, \mathrm{d}\xi_0}{\int \mathcal{V}' f_\mathrm{RE} \, \mathrm{d}\xi_0} 
   = \frac{ \int \mathcal{V}' \{A_p\} e^{-A(p)g(\xi_0)} \, \mathrm{d}\xi_0}{\int \mathcal{V}' e^{-A(p)g(\xi_0)} \, \mathrm{d}\xi_0}.

In the DREAM calculation we account for the contributions from electric field acceleration, synchrotron losses
and slowing down (collisional with screening effect plus mean-force bremsstrahlung loss) using the exact 
advection terms that would be included in a regular ``SUPERTHERMAL`` simulation. 

The critical effective field :math:`E_c^\mathrm{eff}` is then defined as the minimum value of the electric field 
for which there exists a real solution to :math:`U(p) = 0`. That is,

.. math::

   E_c^\mathrm{eff} = \mathrm{min}\Biggl( \frac{\left\langle\boldsymbol{E}\cdot\boldsymbol{B}\right\rangle}{\sqrt{\left\langle B^2 \right\rangle}}  ~ \Bigg| ~ U(p) = 0  \Biggr)

Algorithm
^^^^^^^^^
The critical field calculation has two main steps: initialisation of splines (typically once per simulation), 
and solution of the optimization problem for :math:`E_c^\mathrm{eff}` (at each time step or iteration of the solver).

Initialisation of splines
*************************
Since flux surface averages are significantly more expensive to evaluate than splines, and it is typically sufficient 
for us to resolve the critical field with a relative tolerance of approximately :math:`10^{-3}`, it has proven effective
to create 1D splines (one at each radius) over quantities involving flux-surface or bounce averages (where relatively 
sparse splines are sufficient, sampling only tens of points). 

First, we spline the exponent of the pitch distribution :math:`\xi_0 / \langle \xi \rangle` on a uniform :math:`\xi_0` 
grid. This allows the :math:`g` function to be evaluated efficiently using routines for exact integration of a spline.

Secondly, we identify that all advection terms of interest can be factorised on the form :math:`A_p = a_p(p) \hat{A}_p(\xi_0,\,\theta)` 
where the prefactor depends only on momentum and the remainder only on pitch and poloidal angle. Then, the bounce-average 
of the pitch-dependent part of each advection term :math:`\{\hat{A}_p\}` is splined on a uniform pitch grid in the interval 
:math:`\xi_0 \in [0,1]` since all advection terms of interest are either symmetric or anti-symmetric.

Finally, these steps allow rapid calculation of a spline of the distribution-bounce average of :math:`\hat{A}_p` for each
advection term, sampled in the mapped variable :math:`X = A^2/(1+A)^2 \in [0,\,1]` (in which the functions are 
smoothly varying all the way up to the limit of :math:`A=\infty` corresponding to all runaways having :math:`\xi=1`), 
where :math:`A` is the pitch distribution width parameter appearing in the exponent of :math:`f_\mathrm{RE}`.

This procedure allows the geometric part of the acceleration function :math:`U(p)` to be evaluated essentially for free,
with only the momentum (and unknown-quantity) dependent prefactors requiring further computation, most of which is spent 
on the slowing-down frequency :math:`\nu_s`.


Optimization problem
********************
We solve for the critical field as a nested optimization problem, where the outer layer is the 
one-dimensional root finding problem 

.. math::
   U_e\left(\left\langle \boldsymbol{E}\cdot\boldsymbol{B} \right\rangle\right) = 0,

where :math:`U_e` is the extremum of :math:`U(p)` with respect to :math:`p` 
at a given :math:`\left\langle \boldsymbol{E}\cdot\boldsymbol{B} \right\rangle`. This is solved 
with an unbounded secant method (Newton's method with the numerical derivative taken between 
successive iterations), where the initial guess in each solve is that :math:`E_c^\mathrm{eff}/E_c^\mathrm{tot}`
is constant in time.

The inner layer is the minimization problem 

.. math::
   U_e = \mathrm{min}_p (-U(p))

where the electric field is provided from the outer layer. This is solved using the bounded Brent's method
implemented in the ``gsl_root`` library, requiring an interval in :math:`p` to be given which contains the minimum.
This is done by choosing a narrow interval :math:`p_\mathrm{opt} \times (1 \pm 0.02)`, where :math:`p_\mathrm{opt}`
is the optimum from the previous solve. If the interval does not contain the minimum, it is expanded until a minimum 
is enclosed -- this is typically needed less than once in a thousand solves.


Settings
^^^^^^^^
The model used for :math:`E_c^\mathrm{eff}` is controlled with the following settings:

+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| Option                             | Description                                                                                                                                       |
+====================================+===================================================================================================================================================+
| ``COLLQTY_ECEFF_MODE_EC_TOT``      | Use the approximation :math:`E_c^\mathrm{eff} = E_c^\mathrm{tot}`                                                                                 |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``COLLQTY_ECEFF_MODE_CYLINDRICAL`` | Uses the analytical model Eq (23-24) in `Hesslow PPCF (2018) <https://doi.org/10.1088/1361-6587/aac33e>`_                                         |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``COLLQTY_ECEFF_MODE_SIMPLE``      | ``FULL`` but replacing :math:`\frac{\xi_0'}{\langle \xi'\rangle} = 1` for passing and :math:`0` for trapped electrons in the pitch distribution.  |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``COLLQTY_ECEFF_MODE_FULL``        | Full model outlined above                                                                                                                         |
+------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

.. note::
   In typical scenarios with toroidal geometry, we observe discrepencies 
   of up to 2% between modes ``FULL`` and ``SIMPLE``, less than 10% 
   with ``CYLINDRICAL`` and over 50% with ``EC_TOT``.

.. note::
   ``SIMPLE`` takes essentially the same computation time as ``FULL``,
   and is therefore not recommended except for benchmarking. Compared 
   with the simpler models ``EC_TOT`` and ``CYLINDRICAL``, ``FULL`` can
   make especially fluid simulations substantially slower.

Example
^^^^^^^
An example of how the mode for the critical effective field can be set to ``CYLINDRICAL`` is given below:


.. code-block:: python 

   import DREAM.Settings.Equations.RunawayElectons as Runaways

   ds = DREAMSettings()

   ds.eqsys.n_re.setEceff(RunawayElectrons.COLLQTY_ECEFF_MODE_CYLINDRICAL)




Class documentation
-------------------
.. autoclass:: DREAM.Settings.Equations.RunawayElectrons.RunawayElectrons
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__


