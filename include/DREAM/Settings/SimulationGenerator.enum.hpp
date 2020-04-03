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
    UQTY_N_COLD_EQN_PRESCRIBED=1
};

/////////////////////////////////////
///
/// COLLISION QUANTITY HANDLER SETTINGS
///
/////////////////////////////////////
enum collqty_nu_s_type {
    COLLQTY_NU_S_TYPE_NON_SCREENED=1
};

