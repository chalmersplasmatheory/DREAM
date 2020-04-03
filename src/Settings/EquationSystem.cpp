
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

    // Define unknowns
    DefineUnknowns(eqsys, s, fluidGrid, hottailGrid, runawayGrid);

    // Define equations
    DefineEquations(eqsys, s, fluidGrid, hottailGrid, runawayGrid);

    return eqsys;
}

/**
 * Set the equations of the equation system.
 *
 * eqsys:       Equation system to define quantities in.
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
void SimulationGenerator::DefineEquations(
    EquationSystem *eqsys, Settings *s, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) {
    
}

/**
 * Define the unknowns of the equation system.
 *
 * eqsys:       Equation system to define quantities in.
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
void SimulationGenerator::DefineUnknowns(
    EquationSystem *eqsys, Settings *s, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) {
    // Fluid quantities
    eqsys->SetUnknown(UQTY_E_FIELD, EquationSystem::REGION_FLUID);
    eqsys->SetUnknown(UQTY_N_COLD, fluidGrid);

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        eqsys->SetUnknown(UQTY_F_HOT, hottailGrid);
    } else {
        eqsys->SetUnknown(UQTY_N_HOT, fluidGrid);
    }

    // Runaway quantities
    if (runawayGrid != nullptr) {
        eqsys->SetUnknown(UQTY_F_RE, runawayGrid);
    } else {
        eqsys->SetUnknown(UQTY_N_RE, fluidGrid);
    }
}

