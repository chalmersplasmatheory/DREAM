/**
 * Definition of OptionConstants constants.
 */

#include "DREAM/Settings/OptionConstants.hpp"

using namespace DREAM;

/**
 * NAMES OF UNKNOWN QUANTITIES
 */
const char *OptionConstants::UQTY_E_FIELD           = "E_field";
const char *OptionConstants::UQTY_F_HOT             = "f_hot";
const char *OptionConstants::UQTY_F_RE              = "f_re";
const char *OptionConstants::UQTY_ION_SPECIES       = "n_i";
const char *OptionConstants::UQTY_I_P               = "I_p";
const char *OptionConstants::UQTY_I_WALL            = "I_wall";
const char *OptionConstants::UQTY_J_HOT             = "j_hot";
const char *OptionConstants::UQTY_J_OHM             = "j_ohm";
const char *OptionConstants::UQTY_J_RE              = "j_re";
const char *OptionConstants::UQTY_J_TOT             = "j_tot";
const char *OptionConstants::UQTY_N_COLD            = "n_cold";
const char *OptionConstants::UQTY_N_HOT             = "n_hot";
const char *OptionConstants::UQTY_N_RE              = "n_re";
const char *OptionConstants::UQTY_N_TOT             = "n_tot";
const char *OptionConstants::UQTY_POL_FLUX          = "psi_p";
const char *OptionConstants::UQTY_PSI_WALL          = "psi_wall";
const char *OptionConstants::UQTY_PSI_EDGE          = "psi_edge";
const char *OptionConstants::UQTY_T_COLD            = "T_cold";
const char *OptionConstants::UQTY_V_LOOP_WALL       = "V_loop_w";
const char *OptionConstants::UQTY_W_COLD            = "W_cold";

/**
 * DESCRIPTIONS OF UNKNOWN QUANTITIES
 */
const char *OptionConstants::UQTY_E_FIELD_DESC      = "Parallel electric field <E*B>/sqrt(<B^2>) [V/m]";
const char *OptionConstants::UQTY_F_HOT_DESC        = "Hot electron distribution function [m^-3]";
const char *OptionConstants::UQTY_F_RE_DESC         = "Runaway electron distribution function [m^-3]";
const char *OptionConstants::UQTY_ION_SPECIES_DESC  = "Ion density [m^-3]";
const char *OptionConstants::UQTY_I_P_DESC          = "Total toroidal plasma current [A]";
const char *OptionConstants::UQTY_I_WALL_DESC       = "Wall current [A]";
const char *OptionConstants::UQTY_J_HOT_DESC        = "Hot electron parallel current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_J_OHM_DESC        = "Ohmic current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_J_RE_DESC         = "Runaway electron current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_J_TOT_DESC        = "Total current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_N_COLD_DESC       = "Cold electron density [m^-3]";
const char *OptionConstants::UQTY_N_HOT_DESC        = "Hot electron density [m^-3]";
const char *OptionConstants::UQTY_N_RE_DESC         = "Runaway electron density <n> [m^-3]";
const char *OptionConstants::UQTY_N_TOT_DESC        = "Total electron density [m^-3]";
const char *OptionConstants::UQTY_POL_FLUX_DESC     = "Poloidal magnetic flux normalized to major radius R0 [Vs/m]";
const char *OptionConstants::UQTY_PSI_WALL_DESC     = "Poloidal magnetic flux on tokamak wall (r=rwall) normalized to major radius R0 [Vs/m]";
const char *OptionConstants::UQTY_PSI_EDGE_DESC     = "Poloidal magnetic flux at plasma edge (r=rmax) normalized to major radius R0 [Vs/m]";
const char *OptionConstants::UQTY_T_COLD_DESC       = "Cold electron temperature [eV]";
const char *OptionConstants::UQTY_V_LOOP_WALL_DESC  = "Loop voltage on tokamak wall normalized to major radius R0 [V/m]";
const char *OptionConstants::UQTY_W_COLD_DESC       = "Cold electron total energy density (kinetic plus binding) [J/m^3]";
