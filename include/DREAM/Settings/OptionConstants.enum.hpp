/**
 * This file contains definitions of the enums that are used
 * to indicate options in DREAM. The file is included into the
 * definition of the 'SimulationGenerator' class.
 */

/////////////////////////////////////
///
/// INPUT DATA OPTIONS
///
/////////////////////////////////////
// Interpolation using our own 'Interpolator1D' object
enum prescribed_data_interp {
    // We start from 0 here to remain somewhat compatible
    // with the GSL interpolation interface
    PRESCRIBED_DATA_INTERP_NEAREST=0,
    PRESCRIBED_DATA_INTERP_LINEAR=1,
    PRESCRIBED_DATA_INTERP_LOGARITHMIC=2
};
enum prescribed_data_interp_gsl {
    PRESCRIBED_DATA_INTERP_GSL_LINEAR=1,
    PRESCRIBED_DATA_INTERP_GSL_POLYNOMIAL=2,
    PRESCRIBED_DATA_INTERP_GSL_CSPLINE=3,
    PRESCRIBED_DATA_INTERP_GSL_AKIMA=4,
    PRESCRIBED_DATA_INTERP_GSL_STEFFEN=5
};
enum prescribed_data_interp3d {
    PRESCRIBED_DATA_INTERP3D_NEAREST=0,
    PRESCRIBED_DATA_INTERP3D_LINEAR=1,
    PRESCRIBED_DATA_INTERP3D_LOGARITHMIC=2
};

enum ion_data_type {
    ION_DATA_PRESCRIBED=1,
    ION_DATA_EQUILIBRIUM=2,
    ION_DATA_TYPE_DYNAMIC=3
};

enum ion_opacity_mode {
	OPACITY_MODE_TRANSPARENT=1,
	OPACITY_MODE_GROUND_STATE_OPAQUE=2
};

enum ion_charged_diffusion_mode {
	ION_CHARGED_DIFFUSION_MODE_NONE=1,
	ION_CHARGED_DIFFUSION_MODE_PRESCRIBED=2
};

enum ion_neutral_diffusion_mode {
	ION_NEUTRAL_DIFFUSION_MODE_NONE=1,
	ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED=2
};

enum ion_charged_advection_mode {
	ION_CHARGED_ADVECTION_MODE_NONE=1,
	ION_CHARGED_ADVECTION_MODE_PRESCRIBED=2
};

enum ion_neutral_advection_mode {
	ION_NEUTRAL_ADVECTION_MODE_NONE=1,
	ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED=2
};

enum ion_source_type {
	ION_SOURCE_NONE=1,
	ION_SOURCE_PRESCRIBED=2
};

// Interpolation method for ADAS rate coefficients
enum adas_interp_type {
    ADAS_INTERP_BILINEAR=1,
    ADAS_INTERP_BICUBIC=2
};

/////////////////////////////////////
///
/// GRID OPTIONS
///
/////////////////////////////////////
// Type of radial grid
enum radialgrid_type {
    RADIALGRID_TYPE_CYLINDRICAL=1,
    RADIALGRID_TYPE_TOROIDAL_ANALYTICAL=2,
    RADIALGRID_TYPE_NUMERICAL=3
};

// Type of momentum grid
enum momentumgrid_type {
    MOMENTUMGRID_TYPE_PXI=1,
    MOMENTUMGRID_TYPE_PPARPPERP=2
};

// Type of p grid
enum pxigrid_ptype {
    PXIGRID_PTYPE_UNIFORM=1,
    PXIGRID_PTYPE_BIUNIFORM=2,
    PXIGRID_PTYPE_CUSTOM=3
};

// Type of xi grid
enum pxigrid_xitype {
    PXIGRID_XITYPE_UNIFORM=1,
    PXIGRID_XITYPE_BIUNIFORM=2,
    PXIGRID_XITYPE_UNIFORM_THETA=3,
    PXIGRID_XITYPE_BIUNIFORM_THETA=4,
    PXIGRID_XITYPE_CUSTOM=5,
    PXIGRID_XITYPE_TRAPPED=6
};

// Type of advection interpolation coefficient for jacobian
enum adv_jacobian_mode {
    AD_INTERP_JACOBIAN_LINEAR=1, // does not include non-linear jacobian from flux limiter
    AD_INTERP_JACOBIAN_FULL=2,   // includes non-linear jacobian from flux limiter
    AD_INTERP_JACOBIAN_UPWIND=3  // uses upwind interpolation in the jacobian
};

enum radialgrid_numeric_format {
    RADIALGRID_NUMERIC_FORMAT_LUKE=1
};

/////////////////////////////////////
///
/// SOLVER OPTIONS
///
/////////////////////////////////////
enum solver_type {
    SOLVER_TYPE_LINEARLY_IMPLICIT=1,
    SOLVER_TYPE_NONLINEAR=2
};
// Linear solver type (used by both the linear-implicit
// and nonlinear solvers)
enum linear_solver {
    LINEAR_SOLVER_NONE=0,       // only for backup solver
    LINEAR_SOLVER_LU=1,
    LINEAR_SOLVER_MUMPS=2,
    LINEAR_SOLVER_MKL=3,
    LINEAR_SOLVER_SUPERLU=4,
    LINEAR_SOLVER_GMRES=5
};

/////////////////////////////////////
///
/// TIME STEPPER OPTIONS
///
/////////////////////////////////////
enum timestepper_type {
    TIMESTEPPER_TYPE_CONSTANT=1,
    TIMESTEPPER_TYPE_ADAPTIVE=2,
	TIMESTEPPER_TYPE_IONIZATION=3,
    TIMESTEPPER_TYPE_PYTHON_TERMINATE=4
};

/////////////////////////////////////
///
/// CONDUCTIVITY OPTIONS
///
/////////////////////////////////////
enum conductivity_mode {
    CONDUCTIVITY_MODE_BRAAMS = 1,
    CONDUCTIVITY_MODE_SAUTER_COLLISIONLESS = 2,
    CONDUCTIVITY_MODE_SAUTER_COLLISIONAL = 3
};

enum corrected_conductivity {
    CORRECTED_CONDUCTIVITY_DISABLED = 1,
    CORRECTED_CONDUCTIVITY_ENABLED = 2
};

/////////////////////////////////////
///
/// CURRENT DENSITY OPTIONS
///
/////////////////////////////////////
enum current_profile_type {
	CURRENT_PROFILE_TYPE_J_PARALLEL = 1,
	CURRENT_PROFILE_TYPE_J_DOT_GRADPHI = 2,
	CURRENT_PROFILE_TYPE_JTOR_OVER_R = 3
};

/////////////////////////////////////
///
/// UNKNOWN QUANTITY OPTIONS
///
/////////////////////////////////////
enum uqty_E_field_eqn {
    UQTY_E_FIELD_EQN_PRESCRIBED=1,     // E_field is prescribed by the user
    UQTY_E_FIELD_EQN_SELFCONSISTENT=2, // E_field is prescribed by the user
	UQTY_E_FIELD_EQN_PRESCRIBED_CURRENT=3,// j_ohm is prescribed by the user
};

enum uqty_f_re_inittype {
    UQTY_F_RE_INIT_FORWARD=1,           // Put all particles in p=pMin, xi=+/-1 (sign depending on E)
    UQTY_F_RE_INIT_XI_NEGATIVE=2,       // Put all particles in p=pMin, xi=-1
    UQTY_F_RE_INIT_XI_POSITIVE=3,       // Put all particles in p=pMin, xi=+1
    UQTY_F_RE_INIT_ISOTROPIC=4,         // Distribute all particles isotropically in p=pMin
    UQTY_F_RE_INIT_AVALANCHE=5,			// Distribute particles according to an analytical avalanche distribution
    UQTY_F_RE_INIT_PRESCRIBED=6			// Set prescribed distribution through f_re.setInitialValue
};

enum uqty_V_loop_wall_eqn {
    UQTY_V_LOOP_WALL_EQN_PRESCRIBED=1,     // V_loop on wall (r=b) is prescribed by the user
    UQTY_V_LOOP_WALL_EQN_SELFCONSISTENT=2, // V_loop on wall is evolved self-consistently
    UQTY_V_LOOP_WALL_EQN_TRANSFORMER=3     // V_loop on wall is evolved self-consistently AND externally applied via transformer action
};

enum uqty_n_cold_eqn {
    UQTY_N_COLD_EQN_PRESCRIBED=1,    // n_cold is calculated from ion species as sum_i n_i Z_i
    UQTY_N_COLD_EQN_SELFCONSISTENT=2 // n_cold is calculated from charge neutrality as sum_i n_i Z_i  - n_hot - n_RE
};

enum uqty_T_cold_eqn {
    UQTY_T_COLD_EQN_PRESCRIBED=1,   // T_cold prescribed by the user
    UQTY_T_COLD_SELF_CONSISTENT=2   // T_cold calculated self-consistently
};

enum uqty_T_abl_eqn {
    UQTY_T_ABL_EQN_PRESCRIBED=1,   // T_abl prescribed by the user
    UQTY_T_ABL_SELF_CONSISTENT=2   // T_abl calculated self-consistently
};

enum uqty_T_i_eqn {
    UQTY_T_I_NEGLECT=1,             // Ion temperature not modelled
    UQTY_T_I_INCLUDE=2              // Ion temperature(s) calculated self-consistently
};

enum uqty_distribution_mode {
    UQTY_DISTRIBUTION_MODE_NUMERICAL=1,    // distribution modelled numerically on a kinetic grid
    UQTY_DISTRIBUTION_MODE_ANALYTICAL=2,   // distribution modelled with analytical distribution function
	UQTY_DISTRIBUTION_MODE_PRESCRIBED=3    // distribution is prescribed in time from user input
};

enum uqty_f_hot_dist_mode {                     // Model used for analytic hottail distribution
    UQTY_F_HOT_DIST_MODE_NONREL = 1             // Smith & Verwichte (2008) equation (9-10)
};


/////////////////////////////////////
///
/// COLLISION QUANTITY HANDLER SETTINGS
///
/////////////////////////////////////
enum collqty_lnLambda_type {             // The Coulomb logarithm is...
    COLLQTY_LNLAMBDA_CONSTANT=1,         // the relativistic lnLambda, lnL = lnLc
    COLLQTY_LNLAMBDA_ENERGY_DEPENDENT=2, // energy dependent, separate for collisions with electrons and ions
    COLLQTY_LNLAMBDA_THERMAL=3,          // the thermal lnLambda, lnL = lnLT
    COLLQTY_LNLAMBDA_ION_ION=4           // the ion-ion lnLambda, lnL = lnLii
};

enum collqty_collfreq_mode {
    COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL=1,      // Collision frequencies given in limit T->0 (except in Coulomb logarithm where T_cold enters)
    COLLQTY_COLLISION_FREQUENCY_MODE_FULL=2,              // Full collision frequencies (with Chandrasekhar and error functions etc in the non-relativistic case)
    COLLQTY_COLLISION_FREQUENCY_MODE_ULTRA_RELATIVISTIC=3 // Collision frequencies given in the limit p>>mc, i.e. ignoring the 1/v^2 behavior but keeping the logarithmic increase with gamma
};

enum collqty_collfreq_type {
    COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED=1,           // only free electrons contribute
    COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED=2,                  // free and bound electrons contribute equally
    COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED=3,            // bound electrons contribute via mean excitation energies etc
    COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED_WALKOWIAK=4   // bound electrons contribution with Walkowiak model https://doi.org/10.1063/5.0075859
};

enum collqty_pstar_mode {                // Runaway growth rates are determined from dynamics that are
    COLLQTY_PSTAR_MODE_COLLISIONAL = 1,  // collisional (no trapping correction)
    COLLQTY_PSTAR_MODE_COLLISIONLESS = 2 // collisionless (with trapping correction)
};

enum collqty_screened_diffusion_mode {              // The energy diffusion frequency due to bound electrons are
    COLLQTY_SCREENED_DIFFUSION_MODE_ZERO = 1,       // set to zero
    COLLQTY_SCREENED_DIFFUSION_MODE_MAXWELLIAN = 2  // such that equilibrium distribution is Maxwellian
};

enum collqty_Eceff_mode {
    COLLQTY_ECEFF_MODE_EC_TOT = 1,      // Gives Ectot including all bound electrons (or Ec_free if no impurities/complete screening)
    COLLQTY_ECEFF_MODE_CYLINDRICAL = 2, // Sets Eceff using the Hesslow formula ignoring trapping effects.
    COLLQTY_ECEFF_MODE_SIMPLE = 3,      // An approximate numerical calculation with a simplified account of trapping effects
    COLLQTY_ECEFF_MODE_FULL = 4         // Full 'Lehtinen theory' expression.
};




/////////////////////////////////////
///
/// EQUATION TERM OPTIONS
///
/////////////////////////////////////

enum eqterm_avalanche_mode {                        // Avalanche generation is...
    EQTERM_AVALANCHE_MODE_NEGLECT = 1,              // neglected
    EQTERM_AVALANCHE_MODE_FLUID = 2,                // modelled with fluid growth rate formula
    EQTERM_AVALANCHE_MODE_FLUID_HESSLOW = 3,        // modelled with fluid growth rate formula published by Hesslow et al NF 2019
    EQTERM_AVALANCHE_MODE_KINETIC = 4               // modelled kinetically with RP avalanche source
};

enum eqterm_nonlinear_mode {                        // Non-linear self-collisions are...
    EQTERM_NONLINEAR_MODE_NEGLECT = 1,              // neglected
    EQTERM_NONLINEAR_MODE_NON_REL_ISOTROPIC = 2,    // accounted for with isotropic Landau-Fokker-Planck operator
    EQTERM_NONLINEAR_MODE_NORSEPP = 3               // included with full NORSE++ formalism
};

enum eqterm_bremsstrahlung_mode {                   // Bremsstrahlung radiation reaction is...
    EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT=1,           // neglected
    EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER=2,    // accounted for with an effective force F_br(p)
    EQTERM_BREMSSTRAHLUNG_MODE_BOLTZMANN=3          // accounted for with a linear (Boltzmann) integral operator
};

enum eqterm_ripple_mode {                           // Magnetic ripple pitch scattering
    EQTERM_RIPPLE_MODE_NEGLECT=1,                   // neglected
    EQTERM_RIPPLE_MODE_BOX=2,                       // included with sharp resonance width
    EQTERM_RIPPLE_MODE_GAUSSIAN=3,                  // included with gaussian resonance region
};

enum eqterm_synchrotron_mode {                      // Synchrotron radiation reaction is...
    EQTERM_SYNCHROTRON_MODE_NEGLECT=1,              // neglected
    EQTERM_SYNCHROTRON_MODE_INCLUDE=2               // included
};

enum eqterm_timevaryingb_mode {						// Pitch angle advection due to time-varying B...
	EQTERM_TIMEVARYINGB_MODE_NEGLECT=1,				// neglected
	EQTERM_TIMEVARYINGB_MODE_INCLUDE=2				// included
};

enum eqterm_dreicer_mode {
    EQTERM_DREICER_MODE_NONE=1,                     // Disable Dreicer generation
    EQTERM_DREICER_MODE_CONNOR_HASTIE_NOCORR=2,     // Dreicer based on Connor-Hastie formula (without corrections)
    EQTERM_DREICER_MODE_CONNOR_HASTIE=3,            // Dreicer based on Connor-Hastie formula
    EQTERM_DREICER_MODE_NEURAL_NETWORK=4            // Dreicer using neural network by Hesslow et al
};

enum eqterm_compton_mode {
    EQTERM_COMPTON_MODE_NEGLECT=1,                  // No Compton source
    EQTERM_COMPTON_MODE_FLUID=2,                    // Fluid Compton generation rate
    EQTERM_COMPTON_MODE_KINETIC=3,                  // Kinetic Compton source
};

enum eqterm_frozen_current_mode {
	EQTERM_FROZEN_CURRENT_MODE_DISABLED=1,			// Disable the frozen current mode transport
	EQTERM_FROZEN_CURRENT_MODE_CONSTANT=2,			// Assume momentum-independent radial transport
	EQTERM_FROZEN_CURRENT_MODE_BETAPAR=3			// Assume v_|| scaling of radial transport
};

enum eqterm_transport_bc {
    EQTERM_TRANSPORT_BC_CONSERVATIVE=1,             // Conservative boundary condition at r=rmax (no particles can leave the plasma)
    EQTERM_TRANSPORT_BC_F_0=2,                      // Enforce f = 0 at r > rmax
    EQTERM_TRANSPORT_BC_DF_CONST=3,                  // Assume d^2 f / dr^2 = 0 at r > rmax
};

enum eqterm_ionization_mode {                       // Ionization is modelled with...
    EQTERM_IONIZATION_MODE_FLUID=1,                 // fluid ADAS rate coefficients
    EQTERM_IONIZATION_MODE_KINETIC=2,               // kinetic model
    EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC=3,    // kinetic model with approximate jacobian
    EQTERM_IONIZATION_MODE_FLUID_RE=4,              // approximation of the kinetic ionization rate assuming a mono-energetic RE distribution at 20mc, can be used in both fluid and kinetic mode
};

enum eqterm_particle_source_mode {                  // Equation used for S_particle (the kinetic particle source)
    EQTERM_PARTICLE_SOURCE_ZERO     = 1,            // S_particle = 0
    EQTERM_PARTICLE_SOURCE_IMPLICIT = 2,            // S_particle determined implicitly from density conservation
    EQTERM_PARTICLE_SOURCE_EXPLICIT = 3             // S_particle set explicitly as sum of equation terms that alter electron density
};

enum eqterm_spi_velocity_mode {
    EQTERM_SPI_VELOCITY_MODE_NONE=1,
    EQTERM_SPI_VELOCITY_MODE_PRESCRIBED=2
};

enum eqterm_spi_ablation_mode {
    EQTERM_SPI_ABLATION_MODE_NEGLECT=1,
    EQTERM_SPI_ABLATION_MODE_FLUID_NGS=2,
    EQTERM_SPI_ABLATION_MODE_KINETIC_NGS=3,
    EQTERM_SPI_ABLATION_MODE_NGPS=4
};

enum eqterm_spi_deposition_mode {
    EQTERM_SPI_DEPOSITION_MODE_NEGLECT=1,
    EQTERM_SPI_DEPOSITION_MODE_LOCAL=2,
    EQTERM_SPI_DEPOSITION_MODE_LOCAL_LAST_FLUX_TUBE=3,
    EQTERM_SPI_DEPOSITION_MODE_LOCAL_GAUSSIAN=4
};

enum eqterm_spi_shift_mode {
    EQTERM_SPI_SHIFT_MODE_NEGLECT=1,
    EQTERM_SPI_SHIFT_MODE_PRESCRIBED=2,
    EQTERM_SPI_SHIFT_MODE_ANALYTICAL=3
};

enum eqterm_spi_heat_absorbtion_mode {
    EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT=1,
    EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS=2,
    EQTERM_SPI_HEAT_ABSORBTION_MODE_LOCAL_FLUID_NGS_GAUSSIAN=3
};

enum eqterm_spi_cloud_radius_mode {
    EQTERM_SPI_CLOUD_RADIUS_MODE_NEGLECT=1,
    EQTERM_SPI_CLOUD_RADIUS_MODE_PRESCRIBED_CONSTANT=2,
    EQTERM_SPI_CLOUD_RADIUS_MODE_SELFCONSISTENT=3
};

enum eqterm_spi_abl_ioniz_mode {
    EQTERM_SPI_ABL_IONIZ_MODE_NEUTRAL=1,
    EQTERM_SPI_ABL_IONIZ_MODE_SINGLY_IONIZED=2,
    EQTERM_SPI_ABL_IONIZ_MODE_SELF_CONSISTENT=3
};

enum eqterm_spi_magnetic_field_dependence_mode {
    EQTERM_SPI_MAGNETIC_FIELD_DEPENDENCE_MODE_NEGLECT=1,
    EQTERM_SPI_MAGNETIC_FIELD_DEPENDENCE_MODE_JOREK=2,
};

enum eqterm_particle_source_shape {
    EQTERM_PARTICLE_SOURCE_SHAPE_MAXWELLIAN = 1,    // Maxwellian shape with temperature T_cold
    EQTERM_PARTICLE_SOURCE_SHAPE_DELTA = 2          // Delta function in p=0
};

enum eqterm_hottail_mode {                          // Mode used for hottail runaway generation
    EQTERM_HOTTAIL_MODE_DISABLED = 1,               // Hottail RE generation neglected
    EQTERM_HOTTAIL_MODE_ANALYTIC = 2,               // Ida's MSc thesis (4.24), roughly equivalent to Smith & Verwicthe 2008 Eq (4)
    EQTERM_HOTTAIL_MODE_ANALYTIC_ALT_PC = 3,        // Ida's MSc thesis (4.39)
};

enum eqterm_lcfs_loss_mode {                        // Loss term
    EQTERM_LCFS_LOSS_MODE_DISABLED = 1,
    EQTERM_LCFS_LOSS_MODE_FLUID = 2,
    EQTERM_LCFS_LOSS_MODE_KINETIC = 3};

enum eqterm_tritium_mode {                        // Tritium generation is...
    EQTERM_TRITIUM_MODE_NEGLECT = 1,              // neglected
    EQTERM_TRITIUM_MODE_FLUID = 2,                // Fluid tritium generation rate
    EQTERM_TRITIUM_MODE_KINETIC = 3               // Kinetic tritium generation rate
};


// Option for which parameter to do the 1D interpolation (time or plasma current)
enum svensson_interp1d_param {
    SVENSSON_INTERP1D_TIME=1,
    SVENSSON_INTERP1D_IP=2
};
