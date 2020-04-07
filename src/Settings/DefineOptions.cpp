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
    DefineOptions_EquationSystem(s);
    DefineOptions_RadialGrid(s);
    DefineOptions_HotTailGrid(s);
    DefineOptions_RunawayGrid(s);
    DefineOptions_Solver(s);
    DefineOptions_TimeStepper(s);
}

