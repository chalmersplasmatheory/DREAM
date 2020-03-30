
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;

/**
 * Define the options which can be set for things
 * related to the equation system.
 *
 * s: Settings object to define the options for.
 */
void SimulationGenerator::DefineOptions_EquationSystem(Settings *s) {
    
}

/**
 * Construct an equation system, based on the specification
 * in the given 'Settings' object.
 *
 * s:           Settings object specifying how to construct
 *              the equation system.
 * fluidGrid:   Radial grid for the computation.
 * hottailGrid: Grid on which the hot-tail electron population
 *              is computed.
 * runawayGrid: Grid on which the runaway electron population
 *              is computed.
 *
 * NOTE: The 'hottailGrid' and 'runawayGrid' will be 'nullptr'
 *       if disabled.
 */
EquationSystem *SimulationGenerator::ConstructEquationSystem(
    Settings *s, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) {
    EquationSystem *eqsys = new EquationSystem(fluidGrid, hottailGrid, runawayGrid);

    ////////////////////////////////
    /// DEFINE UNKNOWNS
    ////////////////////////////////
    
    // Fluid quantities
    //eqsys->SetUnknown("E", EquationSystem::REGION_FLUID);
    eqsys->SetUnknown("n_cold", fluidGrid);

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        eqsys->SetUnknown("f_fast", hottailGrid);
    } else {
        eqsys->SetUnknown("n_fast", fluidGrid);
    }

    // Runaway quantities
    if (runawayGrid != nullptr) {
        eqsys->SetUnknown("f_re", runawayGrid);
    } else {
        eqsys->SetUnknown("n_RE", fluidGrid);
    }

    return eqsys;
}

