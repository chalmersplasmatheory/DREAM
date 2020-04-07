
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;


#define EQUATIONSYSTEM "equationsystem"

/**
 * Define the options which can be set for things
 * related to the equation system.
 *
 * s: Settings object to define the options for.
 */
void SimulationGenerator::DefineOptions_EquationSystem(Settings *s) {
    s->DefineSetting(EQUATIONSYSTEM "/n_cold/type", "Type of equation to use for determining the cold electron density", (int_t)UQTY_N_COLD_EQN_PRESCRIBED);
    s->DefineSetting(EQUATIONSYSTEM "/n_cold/data/n", "Prescribed cold electron density", 0, (real_t*)nullptr);
    s->DefineSetting(EQUATIONSYSTEM "/n_cold/data/r", "Radial grid used for prescribing the cold electron density", 0, (real_t*)nullptr);
    s->DefineSetting(EQUATIONSYSTEM "/n_cold/data/t", "Times corresponding to prescribed cold electron density", 0, (real_t*)nullptr);
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

    // Construct the time stepper
    ConstructTimeStepper(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(eqsys, s, fluidGrid, hottailGrid, runawayGrid);

    // Construct equations according to settings
    ConstructEquations(eqsys, s);

    eqsys->ProcessSystem();
    // Construct solver (must be done after processing equation system,
    // since we need to know which unknowns are "non-trivial",
    // i.e. need to show up in the solver matrices)
    ConstructSolver(eqsys, s);

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
void SimulationGenerator::ConstructEquations(
    EquationSystem *eqsys, Settings *s
) {
    // Fluid equations
    ConstructEquation_n_cold(eqsys, s);
}

/**
 * Construct the unknowns of the equation system.
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
void SimulationGenerator::ConstructUnknowns(
    EquationSystem *eqsys, Settings *s, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid *runawayGrid
) {
    // Fluid quantities
    //eqsys->SetUnknown(UQTY_E_FIELD, fluidGrid);
    eqsys->SetUnknown(UQTY_N_COLD, fluidGrid);
    eqsys->SetUnknown(UQTY_N_HOT, fluidGrid);
    eqsys->SetUnknown(UQTY_N_RE, fluidGrid);

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        eqsys->SetUnknown(UQTY_F_HOT, hottailGrid);
    }

    // Runaway quantities
    if (runawayGrid != nullptr) {
        eqsys->SetUnknown(UQTY_F_RE, runawayGrid);
    }
}

