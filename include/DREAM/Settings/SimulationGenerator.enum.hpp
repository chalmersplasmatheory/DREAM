/**
 * This file contains definitions of the enums that are used
 * to indicate options in DREAM. The file is included into the
 * definition of the 'SimulationGenerator' class.
 */

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
/// UNKNOWN QUANTITY OPTIONS
///
/////////////////////////////////////
enum uqty_n_cold_eqn {
    UQTY_N_COLD_EQN_PRESCRIBED=1,     // n_cold is calcaulted from ion species as sum_i n_i Z_i
    UQTY_N_COLD_EQN_SELFCONSISTENT=2, // n_cold is calculated from charge neutrality as sum_i n_i Z_i  - n_hot - n_RE
};

/////////////////////////////////////
///
/// COLLISION QUANTITY HANDLER SETTINGS
///
/////////////////////////////////////
enum collqty_collfreq_type {
    COLLQTY_COLLISION_FREQUENCY_TYPE_SUPERTHERMAL_COMPLETELY_SCREENED=1, // only free electrons contribute 
    COLLQTY_COLLISION_FREQUENCY_TYPE_SUPERTHERMAL_NON_SCREENED=2,        // free and bound electrons contribute equally
    COLLQTY_COLLISION_FREQUENCY_TYPE_SUPERTHERMAL_PARTIALLY_SCREENED=3   // bound electrons contribute via mean excitation energies  
};

