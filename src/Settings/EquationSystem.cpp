
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

    s->DefineSetting(EQUATIONSYSTEM "/n_cold/type", "Type of equation to use for determining the cold electron density", (int_t)OptionConstants::UQTY_N_COLD_EQN_PRESCRIBED);
    DefineDataRT(EQUATIONSYSTEM "/n_cold", s);
    
    s->DefineSetting(EQUATIONSYSTEM "/T_cold/type", "Type of equation to use for determining the electric field evolution", (int_t)OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED);
    DefineDataRT(EQUATIONSYSTEM "/T_cold", s);
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
    enum OptionConstants::momentumgrid_type re_type, FVM::Grid *runawayGrid,
    ADAS *adas
) {
    EquationSystem *eqsys = new EquationSystem(fluidGrid, ht_type, hottailGrid, re_type, runawayGrid);

    // Construct the time stepper
    ConstructTimeStepper(eqsys, s);

    // Construct unknowns
    ConstructUnknowns(eqsys, s, fluidGrid, hottailGrid, runawayGrid);

    

    // Construct equations according to settings
    ConstructEquations(eqsys, s, adas);

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
    EquationSystem *eqsys, Settings *s, ADAS *adas
) {
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();
    FVM::Grid *runawayGrid = eqsys->GetRunawayGrid();
    enum OptionConstants::momentumgrid_type ht_type = eqsys->GetHotTailGridType();
    enum OptionConstants::momentumgrid_type re_type = eqsys->GetRunawayGridType();

    // Fluid equations
    ConstructEquation_Ions(eqsys, s, adas);
    IonHandler *ionHandler = eqsys->GetIonHandler();
    // Construct collision quantity handlers
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    if (hottailGrid != nullptr) {
        CollisionQuantityHandler *cqh = ConstructCollisionQuantityHandler("hottailgrid", ht_type, hottailGrid, unknowns, ionHandler, s);
        eqsys->SetHotTailCollisionHandler(cqh);
    } else if (runawayGrid != nullptr) {
        CollisionQuantityHandler *cqh = ConstructCollisionQuantityHandler("runawaygrid", re_type, runawayGrid, unknowns, ionHandler, s);
        eqsys->SetRunawayCollisionHandler(cqh);
    }

    ConstructEquation_E_field(eqsys, s);
    ConstructEquation_n_cold(eqsys, s);
    ConstructEquation_n_hot(eqsys, s);
    ConstructEquation_T_cold(eqsys, s);

    // Helper quantities
    ConstructEquation_n_tot(eqsys, s);

    // Hot-tail quantities
    if (eqsys->HasHotTailGrid()) {
        ConstructEquation_f_hot(eqsys, s);
    }

    // NOTE: The runaway number may depend explicitly on the
    // hot-tail equation and must therefore be constructed
    // AFTER the call to 'ConstructEquation_f_hot()'
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
    EquationSystem *eqsys, Settings *s, FVM::Grid *fluidGrid,
    FVM::Grid *hottailGrid, FVM::Grid*
) {
    // Fluid quantities
    eqsys->SetUnknown(OptionConstants::UQTY_E_FIELD, fluidGrid);
    eqsys->SetUnknown(OptionConstants::UQTY_N_HOT, fluidGrid);
    // (note: n_cold should be defined after 'n_hot', since when
    // evaluated self-consistently, 'n_cold' depends on 'n_hot')
    eqsys->SetUnknown(OptionConstants::UQTY_N_COLD, fluidGrid);
    eqsys->SetUnknown(OptionConstants::UQTY_N_RE, fluidGrid);
    eqsys->SetUnknown(OptionConstants::UQTY_T_COLD, fluidGrid);

    // Fluid helper quantities
    eqsys->SetUnknown(OptionConstants::UQTY_N_TOT, fluidGrid);

    len_t nIonChargeStates = GetNumberOfIonChargeStates(s);
    eqsys->SetUnknown(OptionConstants::UQTY_ION_SPECIES, fluidGrid, nIonChargeStates);

    // Hot-tail quantities
    if (hottailGrid != nullptr) {
        eqsys->SetUnknown(OptionConstants::UQTY_F_HOT, hottailGrid);
    }

    // Runaway quantities
    /*if (runawayGrid != nullptr) {
        eqsys->SetUnknown(OptionConstants::UQTY_F_RE, runawayGrid);
    }*/
}
