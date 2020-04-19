
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
    s->DefineSetting(EQUATIONSYSTEM "/E_field/type", "Type of equation to use for determining the electric field evolution", (int_t)OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED);
    DefineDataRT(EQUATIONSYSTEM "/E_field", s);

    DefineDataR2P(EQUATIONSYSTEM "/f_hot", s, "init");

    s->DefineSetting(EQUATIONSYSTEM "/n_cold/type", "Type of equation to use for determining the cold electron density", (int_t)OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED);
    DefineDataRT(EQUATIONSYSTEM "/n_cold", s);
}

/**
 * Construct an equation system, based on the specification
 * in the given 'Settings' object.
 *
 * s:           Settings object specifying how to construct
 *              the equation system.
 * fluidGrid:   Radial grid for the computation.
 * ht_type:     Exact type of the hot-tail momentum grid.
 * hottailGrid: Grid on which the hot-tail electron population
 *              is computed.
 * re_type:     Exact type of the runaway momentum grid.
 * runawayGrid: Grid on which the runaway electron population
 *              is computed.
 *
 * NOTE: The 'hottailGrid' and 'runawayGrid' will be 'nullptr'
 *       if disabled.
 */
EquationSystem *SimulationGenerator::ConstructEquationSystem(
    Settings *s, FVM::Grid *fluidGrid,
    enum OptionConstants::momentumgrid_type ht_type, FVM::Grid *hottailGrid,
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid
) {
    EquationSystem *eqsys = new EquationSystem(fluidGrid, ht_type, hottailGrid, re_type, runawayGrid);

    // Construct the time stepper
    ConstructTimeStepper(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(eqsys, s, fluidGrid, hottailGrid, runawayGrid);

    // Construct collision quantity handlers
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    if (hottailGrid != nullptr) {
        CollisionQuantityHandler *cqh = ConstructCollisionQuantityHandler("hottailgrid", ht_type, hottailGrid, unknowns, s);
        eqsys->SetHotTailCollisionHandler(cqh);
    } else {
        CollisionQuantityHandler *cqh = ConstructCollisionQuantityHandler("runawaygrid", re_type, runawayGrid, unknowns, s);
        eqsys->SetRunawayCollisionHandler(cqh);
    }

    // Construct equations according to settings
    ConstructEquations(eqsys, s);

    // Figure out which unknowns must be part of the matrix,
    // and set initial values for those quantities which don't
    // yet have an initial value.
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
    ConstructEquation_E_field(eqsys, s);
    ConstructEquation_n_cold(eqsys, s);
    ConstructEquation_n_hot(eqsys, s);

    // Hot-tail quantities
    if (eqsys->HasHotTailGrid()) {
        ConstructEquation_f_hot(eqsys, s);
    }

    ConstructEquation_n_re(eqsys, s);
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
    EquationSystem *eqsys, Settings* /*s*/, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid*
) {
    // Fluid quantities
    eqsys->SetUnknown(OptionConstants::UQTY_E_FIELD, fluidGrid);
    eqsys->SetUnknown(OptionConstants::UQTY_N_COLD, fluidGrid);
    eqsys->SetUnknown(OptionConstants::UQTY_N_HOT, fluidGrid);
    eqsys->SetUnknown(OptionConstants::UQTY_N_RE, fluidGrid);

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        eqsys->SetUnknown(OptionConstants::UQTY_F_HOT, hottailGrid);
    }

    // Runaway quantities
    /*if (runawayGrid != nullptr) {
        eqsys->SetUnknown(OptionConstants::UQTY_F_RE, runawayGrid);
    }*/
}

