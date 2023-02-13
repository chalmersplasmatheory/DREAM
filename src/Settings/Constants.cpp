/**
 * Definition of OptionConstants constants.
 */

#include "DREAM/Settings/OptionConstants.hpp"

using namespace DREAM;

#define DEF_UQTY(VARNAME,NAME,DESC) \
		const char *OptionConstants::UQTY_ ## VARNAME = NAME; \
		const char *OptionConstants::UQTY_ ## VARNAME ## _DESC = DESC

/**
 * NAMES OF UNKNOWN QUANTITIES
 */
DEF_UQTY(E_FIELD,         "E_field",      "Parallel electric field <E*B>/sqrt(<B^2>) [V/m]");
DEF_UQTY(F_HOT,           "f_hot",        "Hot electron distribution function [m^-3]");
DEF_UQTY(F_RE,            "f_re",         "Runaway electron distribution function [m^-3]");
DEF_UQTY(ION_SPECIES,     "n_i",          "Ion density [m^-3]");
DEF_UQTY(ION_SPECIES_ABL, "n_i_abl",      "Flux surface averaged density of ablated but not yet equilibrated ions [m^-3]");
DEF_UQTY(I_P,             "I_p",          "Total toroidal plasma current [A]");
DEF_UQTY(I_WALL,          "I_wall",       "Wall current [A]");
DEF_UQTY(J_HOT,           "j_hot",        "Hot electron parallel current density j_||*Bmin/B [A/m^2]");
DEF_UQTY(J_OHM,           "j_ohm",        "Ohmic current density j_||*Bmin/B [A/m^2]");
DEF_UQTY(J_RE,            "j_re",         "Runaway electron current density j_||*Bmin/B [A/m^2]");
DEF_UQTY(J_TOT,           "j_tot",        "Total current density j_||*Bmin/B [A/m^2]");
DEF_UQTY(K_I,             "k_I",          "Transport scaling coefficient in frozen current mode");
DEF_UQTY(N_ABL,           "n_abl",        "Flux surface averaged density of ablated but not yet equilibrated electrons [m^-3]");
DEF_UQTY(N_COLD,          "n_cold",       "Cold electron density [m^-3]");
DEF_UQTY(N_HOT,           "n_hot",        "Hot electron density [m^-3]");
DEF_UQTY(N_RE,            "n_re",         "Runaway electron density <n> [m^-3]");
DEF_UQTY(N_RE_NEG,        "n_re_neg",     "Runaway electron density in negative direction <n> [m^-3]");
DEF_UQTY(N_TOT,           "n_tot",        "Total electron density [m^-3]");
DEF_UQTY(NI_DENS,         "N_i",          "Total ion density of each species [m^-3]");
DEF_UQTY(POL_FLUX,        "psi_p",        "Poloidal magnetic flux normalized to major radius R0 [Vs/m]");
DEF_UQTY(PSI_EDGE,        "psi_edge",     "Poloidal magnetic flux at plasma edge (r=rmax) normalized to major radius R0 [Vs/m]");
DEF_UQTY(Q_HOT,           "q_hot",        "Hot out going heat flux density in all directions [J/(s m^2)]");
DEF_UQTY(PSI_TRANS,       "psi_trans",    "Poloidal magnetic flux at transformer normalized to major radius R0 [Vs/m]");
DEF_UQTY(PSI_WALL,        "psi_wall",     "Poloidal magnetic flux on tokamak wall (r=rwall) normalized to major radius R0 [Vs/m]");
DEF_UQTY(S_PARTICLE,      "S_particle",   "Rate at which cold-electron density is added [m^-3 s^-1]");
DEF_UQTY(T_ABL,           "T_abl",        "Ablated but not equilibrated electron temperature [eV]");
DEF_UQTY(T_COLD,          "T_cold",       "Cold electron temperature [eV]");
DEF_UQTY(TAU_COLL,        "tau_coll",     "Time-integrated relativistic collision frequency for analytic hottail formula");
DEF_UQTY(V_LOOP_TRANS,    "V_loop_trans", "Loop voltage applied at transformer normalized to major radius R0 [V/m]");
DEF_UQTY(V_LOOP_WALL,     "V_loop_w",     "Loop voltage on tokamak wall normalized to major radius R0 [V/m]");
DEF_UQTY(V_P,             "v_p",          "Pellet shard velocities (Cartesian) [m/s]");
DEF_UQTY(W_ABL,           "W_abl",        "Flux surface averaged ablated but not equilibrated electron energy density (3n_abl T_abl/2) [J/m^3]");
DEF_UQTY(W_COLD,          "W_cold",       "Cold electron energy density (3nT/2) [J/m^3]");
DEF_UQTY(W_HOT,           "W_hot",        "Hot electron energy density [J/m^3]");
DEF_UQTY(WI_ENER,         "W_i",          "Total ion energy density (3N_iT_i/2) of each species [J/m^3]");
DEF_UQTY(X_P,             "x_p",          "Pellet shard coordinates (Cartesian) [m]");
DEF_UQTY(Y_P,             "Y_p",          "Pellet shard radii to the power of 5/3 [m^(5/3)]");

