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
    PRESCRIBED_DATA_INTERP_LINEAR=1
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
    PRESCRIBED_DATA_INTERP3D_LINEAR=1
};

enum ion_data_type {
    ION_DATA_PRESCRIBED=1,
    ION_DATA_EQUILIBRIUM=2,
    ION_DATA_TYPE_DYNAMIC=3
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
    RADIALGRID_TYPE_TOROIDAL_ANALYTICAL=2
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
    PXIGRID_XITYPE_CUSTOM=5
};

// Type of advection interpolation coefficient for jacobian
enum adv_jacobian_mode {
    AD_INTERP_JACOBIAN_LINEAR=1, // does not include non-linear jacobian from flux limiter 
    AD_INTERP_JACOBIAN_FULL=2,   // includes non-linear jacobian from flux limiter
    AD_INTERP_JACOBIAN_UPWIND=3  // uses upwind interpolation in the jacobian 
};

/////////////////////////////////////
///
/// SOLVER OPTIONS
///
/////////////////////////////////////
enum solver_type {
    SOLVER_TYPE_LINEARLY_IMPLICIT=1,
    SOLVER_TYPE_NONLINEAR=2,
    SOLVER_TYPE_NONLINEAR_SNES=3
};
// Linear solver type (used by both the linear-implicit
// and nonlinear solvers)
enum linear_solver {
    LINEAR_SOLVER_LU=1,
    LINEAR_SOLVER_MUMPS=2
};

/////////////////////////////////////
///
/// TIME STEPPER OPTIONS
///
/////////////////////////////////////
enum timestepper_type {
    TIMESTEPPER_TYPE_CONSTANT=1,
    TIMESTEPPER_TYPE_ADAPTIVE=2
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
/// UNKNOWN QUANTITY OPTIONS
///
/////////////////////////////////////
enum uqty_E_field_eqn {
    UQTY_E_FIELD_EQN_PRESCRIBED=1,     // E_field is prescribed by the user
    UQTY_E_FIELD_EQN_SELFCONSISTENT=2, // E_field is prescribed by the user
};

enum uqty_V_loop_wall_eqn {
    UQTY_V_LOOP_WALL_EQN_PRESCRIBED=1,     // V_loop on wall (r=b) is prescribed by the user
    UQTY_V_LOOP_WALL_EQN_SELFCONSISTENT=2, // V_loop on wall is evolved self-consistently
};

enum uqty_n_cold_eqn {
    UQTY_N_COLD_EQN_PRESCRIBED=1,    // n_cold is calculated from ion species as sum_i n_i Z_i
    UQTY_N_COLD_EQN_SELFCONSISTENT=2 // n_cold is calculated from charge neutrality as sum_i n_i Z_i  - n_hot - n_RE
};

enum uqty_T_cold_eqn {
    UQTY_T_COLD_EQN_PRESCRIBED=1,   // T_cold prescribed by the user
    UQTY_T_COLD_SELF_CONSISTENT=2   // T_cold calculated self-consistently
};

enum uqty_T_i_eqn {
    UQTY_T_I_NEGLECT=1,             // Ion temperature not modelled
    UQTY_T_I_INCLUDE=2              // Ion temperature(s) calculated self-consistently
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
    COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED=1, // only free electrons contribute 
    COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED=2,        // free and bound electrons contribute equally
    COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED=3   // bound electrons contribute via mean excitation energies etc
};

enum collqty_pstar_mode {                // Runaway growth rates are determined from dynamics that are
    COLLQTY_PSTAR_MODE_COLLISIONAL = 1,  // collisional (no trapping correction)
    COLLQTY_PSTAR_MODE_COLLISIONLESS = 2 // collisionless (with trapping correction)
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

enum eqterm_transport_bc {
    EQTERM_TRANSPORT_BC_CONSERVATIVE=1,             // Conservative boundary condition at r=rmax (no particles can leave the plasma)
    EQTERM_TRANSPORT_BC_F_0=2                       // Enforce f = 0 at r > rmax
};

enum eqterm_ionization_mode {                       // Ionization is modelled with...
    EQTERM_IONIZATION_MODE_FLUID=1,                 // fluid ADAS rate coefficients
    EQTERM_IONIZATION_MODE_KINETIC=2,               // kinetic model
    EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC=3,    // kinetic model with approximate jacobian
};

enum eqterm_particle_source_mode {                  // Equation used for S_particle (the kinetic particle source) 
    EQTERM_PARTICLE_SOURCE_ZERO     = 1,            // S_particle = 0
    EQTERM_PARTICLE_SOURCE_IMPLICIT = 2,            // S_particle determined implicitly from density conservation
    EQTERM_PARTICLE_SOURCE_EXPLICIT = 3             // S_particle set explicitly as sum of equation terms that alter electron density
};
