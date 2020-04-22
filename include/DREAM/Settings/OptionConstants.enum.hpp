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
    ION_DATA_FULLY_IONIZED=2,
    ION_DATA_EQUILIBRIUM=3
};

/////////////////////////////////////
///
/// GRID OPTIONS
///
/////////////////////////////////////
// Type of radial grid
enum radialgrid_type {
    RADIALGRID_TYPE_CYLINDRICAL=1
};

// Type of momentum grid
enum momentumgrid_type {
    MOMENTUMGRID_TYPE_PXI=1,
    MOMENTUMGRID_TYPE_PPARPPERP=2
};

// Type of p grid
enum pxigrid_ptype {
    PXIGRID_PTYPE_UNIFORM=1
};

// Type of xi grid
enum pxigrid_xitype {
    PXIGRID_XITYPE_UNIFORM=1
};

/////////////////////////////////////
///
/// SOLVER OPTIONS
///
/////////////////////////////////////
enum solver_type {
    SOLVER_TYPE_LINEARLY_IMPLICIT=1,
    SOLVER_TYPE_NONLINEAR_SNES=2
};

/////////////////////////////////////
///
/// TIME STEPPER OPTIONS
///
/////////////////////////////////////
enum timestepper_type {
    TIMESTEPPER_TYPE_CONSTANT=1
};

/////////////////////////////////////
///
/// UNKNOWN QUANTITY OPTIONS
///
/////////////////////////////////////
enum uqty_E_field_eqn {
    UQTY_E_FIELD_EQN_PRESCRIBED=1   // E_field is prescribed by the user
};

enum uqty_n_cold_eqn {
    UQTY_N_COLD_EQN_PRESCRIBED=1,    // n_cold is calculated from ion species as sum_i n_i Z_i
    UQTY_N_COLD_EQN_SELFCONSISTENT=2 // n_cold is calculated from charge neutrality as sum_i n_i Z_i  - n_hot - n_RE
};

/////////////////////////////////////
///
/// COLLISION QUANTITY HANDLER SETTINGS
///
/////////////////////////////////////
enum collqty_lnLambda_type {             // The Coulomb logarithm is... 
    COLLQTY_LNLAMBDA_CONSTANT=1,         // the relativistic lnLambda, lnL = lnLc
    COLLQTY_LNLAMBDA_ENERGY_DEPENDENT=2  // energy dependent, separate for collisions with electrons and ions
};

enum collqty_collfreq_mode {
    COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL=1, // Collision frequencies given in limit T->0 (except in Coulomb logarithm where T_cold enters)
    COLLQTY_COLLISION_FREQUENCY_MODE_FULL=2          // Full collision frequencies (with Chandrasekhar and error functions etc in the non-relativistic case)
};

enum collqty_collfreq_type {
    COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED=1, // only free electrons contribute 
    COLLQTY_COLLISION_FREQUENCY_TYPE_NON_SCREENED=2,        // free and bound electrons contribute equally
    COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED=3   // bound electrons contribute via mean excitation energies etc
};


/////////////////////////////////////
///
/// EQUATION TERM OPTIONS
///
/////////////////////////////////////

enum eqterm_bremsstrahlung_mode {                // Bremsstrahlung radiation reaction is...
    EQTERM_BREMSSTRAHLUNG_MODE_NEGLECT=1,        // neglected
    EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER=2, // accounted for with an effective force F_br(p)
    EQTERM_BREMSSTRAHLUNG_MODE_BOLTZMANN=3       // accounted for with a linear (Boltzmann) integral operator
};

enum eqterm_synchrotron_mode {         // Synchrotron radiation reaction is...
    EQTERM_SYNCHROTRON_MODE_NEGLECT=1, // neglected 
    EQTERM_SYNCHROTRON_MODE_INCLUDE=2  // included
};


