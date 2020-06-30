/**
 * This file contains the routine 'DefineOptions', which
 * defines all available options in DREAM.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;

/**
 * Define options for the given settings object.
 *
 * s: Settings object to define available options for.
 */
void SimulationGenerator::DefineOptions(Settings *s) {
    DefineOptions_ADAS(s);
    DefineOptions_CollisionQuantityHandler(s);
    DefineOptions_EquationSystem(s);
    DefineOptions_Initializer(s);
    DefineOptions_RadialGrid(s);
    DefineOptions_HotTailGrid(s);
    DefineOptions_ElectricField(s);
    DefineOptions_T_cold(s);
    DefineOptions_f_hot(s);
    DefineOptions_j_ohm(s);
    DefineOptions_j_tot(s);
    DefineOptions_Ions(s);
    DefineOptions_n_re(s);
    DefineOptions_RunawayGrid(s);
    DefineOptions_Solver(s);
    DefineOptions_TimeStepper(s);
    DefineOptions_OtherQuantities(s);
}

/**
 * Define options for the ADAS database.
 */
void SimulationGenerator::DefineOptions_ADAS(Settings *s) {
    s->DefineSetting("/ADAS_interpolation", "Interpolation method for ADAS rate coefficients", (int_t)OptionConstants::ADAS_INTERP_BICUBIC);
}

