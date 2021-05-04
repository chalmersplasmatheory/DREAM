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
const char *OptionConstants::UQTY_ION_SPECIES_ABL   = "n_i_abl";
const char *OptionConstants::UQTY_I_P               = "I_p";
const char *OptionConstants::UQTY_I_WALL            = "I_wall";
const char *OptionConstants::UQTY_J_HOT             = "j_hot";
const char *OptionConstants::UQTY_J_OHM             = "j_ohm";
const char *OptionConstants::UQTY_J_RE              = "j_re";
const char *OptionConstants::UQTY_J_TOT             = "j_tot";
const char *OptionConstants::UQTY_N_ABL             = "n_abl";
const char *OptionConstants::UQTY_N_COLD            = "n_cold";
const char *OptionConstants::UQTY_N_HOT             = "n_hot";
const char *OptionConstants::UQTY_N_RE              = "n_re";
const char *OptionConstants::UQTY_N_TOT             = "n_tot";
const char *OptionConstants::UQTY_NI_DENS           = "N_i";
const char *OptionConstants::UQTY_POL_FLUX          = "psi_p";
const char *OptionConstants::UQTY_PSI_WALL          = "psi_wall";
const char *OptionConstants::UQTY_PSI_EDGE          = "psi_edge";
const char *OptionConstants::UQTY_Q_HOT             = "q_hot";
const char *OptionConstants::UQTY_S_PARTICLE        = "S_particle";
const char *OptionConstants::UQTY_T_ABL             = "T_abl";
const char *OptionConstants::UQTY_T_COLD            = "T_cold";
const char *OptionConstants::UQTY_TAU_COLL          = "tau_coll";
const char *OptionConstants::UQTY_V_LOOP_WALL       = "V_loop_w";
const char *OptionConstants::UQTY_V_P               = "v_p";
const char *OptionConstants::UQTY_W_ABL             = "W_abl";
const char *OptionConstants::UQTY_W_COLD            = "W_cold";
const char *OptionConstants::UQTY_W_HOT             = "W_hot";
const char *OptionConstants::UQTY_WI_ENER           = "W_i";
const char *OptionConstants::UQTY_X_P               = "x_p";
const char *OptionConstants::UQTY_Y_P               = "Y_p";

/**
 * DESCRIPTIONS OF UNKNOWN QUANTITIES
 */
const char *OptionConstants::UQTY_E_FIELD_DESC      = "Parallel electric field <E*B>/sqrt(<B^2>) [V/m]";
const char *OptionConstants::UQTY_F_HOT_DESC        = "Hot electron distribution function [m^-3]";
const char *OptionConstants::UQTY_F_RE_DESC         = "Runaway electron distribution function [m^-3]";
const char *OptionConstants::UQTY_ION_SPECIES_DESC  = "Ion density [m^-3]";
const char *OptionConstants::UQTY_ION_SPECIES_ABL_DESC = "Flux surface averaged density of ablated but not yet equilibrated ions [m^-3]";
const char *OptionConstants::UQTY_I_P_DESC          = "Total toroidal plasma current [A]";
const char *OptionConstants::UQTY_I_WALL_DESC       = "Wall current [A]";
const char *OptionConstants::UQTY_J_HOT_DESC        = "Hot electron parallel current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_J_OHM_DESC        = "Ohmic current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_J_RE_DESC         = "Runaway electron current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_J_TOT_DESC        = "Total current density j_||*Bmin/B [A/m^2]";
const char *OptionConstants::UQTY_N_ABL_DESC        = "Flux surface averaged density of ablated but not yet equilibrated electrons [m^-3]";
const char *OptionConstants::UQTY_N_COLD_DESC       = "Cold electron density [m^-3]";
const char *OptionConstants::UQTY_N_HOT_DESC        = "Hot electron density [m^-3]";
const char *OptionConstants::UQTY_N_RE_DESC         = "Runaway electron density <n> [m^-3]";
const char *OptionConstants::UQTY_N_TOT_DESC        = "Total electron density [m^-3]";
const char *OptionConstants::UQTY_NI_DENS_DESC      = "Total ion density of each species [m^-3]";
const char *OptionConstants::UQTY_POL_FLUX_DESC     = "Poloidal magnetic flux normalized to major radius R0 [Vs/m]";
const char *OptionConstants::UQTY_PSI_WALL_DESC     = "Poloidal magnetic flux on tokamak wall (r=rwall) normalized to major radius R0 [Vs/m]";
const char *OptionConstants::UQTY_PSI_EDGE_DESC     = "Poloidal magnetic flux at plasma edge (r=rmax) normalized to major radius R0 [Vs/m]";
const char *OptionConstants::UQTY_Q_HOT_DESC        = "Hot out going heat flux density in all directions [J/(s m^2)]";
const char *OptionConstants::UQTY_S_PARTICLE_DESC   = "Rate at which cold-electron density is added [m^-3 s^-1]";
const char *OptionConstants::UQTY_T_ABL_DESC        = "Ablated but not equilibrated electron temperature [eV]";
const char *OptionConstants::UQTY_T_COLD_DESC       = "Cold electron temperature [eV]";
const char *OptionConstants::UQTY_TAU_COLL_DESC     = "Time-integrated relativistic collision frequency for analytic hottail formula";
const char *OptionConstants::UQTY_V_LOOP_WALL_DESC  = "Loop voltage on tokamak wall normalized to major radius R0 [V/m]";
const char *OptionConstants::UQTY_V_P_DESC          = "Pellet shard velocities (Cartesian) [m/s]";
const char *OptionConstants::UQTY_W_ABL_DESC        = "Flux surface averaged ablated but not equilibrated electron energy density (3n_abl T_abl/2) [J/m^3]";
const char *OptionConstants::UQTY_W_COLD_DESC       = "Cold electron energy density (3nT/2) [J/m^3]";
const char *OptionConstants::UQTY_W_HOT_DESC        = "Hot electron energy density [J/m^3]";
const char *OptionConstants::UQTY_WI_ENER_DESC      = "Total ion energy density (3N_iT_i/2) of each species [J/m^3]";
const char *OptionConstants::UQTY_X_P_DESC          = "Pellet shard coordinates (Cartesian) [m]";
const char *OptionConstants::UQTY_Y_P_DESC          = "Pellet shard radii to the power of 5/3 [m^(5/3)]";

