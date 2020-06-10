/**
 * Construct a time stepper object.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "DREAM/Solver/SolverSNES.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "solver"


/**
 * Define options for the time stepper.
 * 
 * s: Settings object to define settings in.
 */
void SimulationGenerator::DefineOptions_Solver(Settings *s) {
    s->DefineSetting(MODULENAME "/type", "Equation system solver type", (int_t)OptionConstants::SOLVER_TYPE_NONLINEAR_SNES);
}

/**
 * Construct a Solver object according to the settings.
 *
 * eqsys: Equation system object to assign Solver object to.
 * s:     Settings object specifying how to construct the Solver object.
 */
void SimulationGenerator::ConstructSolver(EquationSystem *eqsys, Settings *s) {
    enum OptionConstants::solver_type type = (enum OptionConstants::solver_type)s->GetInteger(MODULENAME "/type");

    FVM::UnknownQuantityHandler *u = eqsys->GetUnknownHandler();
    vector<UnknownQuantityEquation*> *eqns = eqsys->GetEquations();

    Solver *solver;
    switch (type) {
        case OptionConstants::SOLVER_TYPE_LINEARLY_IMPLICIT:
            solver = ConstructSolver_linearly_implicit(s, u, eqns);
            break;

        case OptionConstants::SOLVER_TYPE_NONLINEAR_SNES:
            solver = ConstructSolver_nonlinear_snes(s, u, eqns);
            break;

        default:
            throw SettingsException(
                "Unrecognized solver type: %d.", type
            );
    }

    eqsys->SetSolver(solver);
    solver->SetCollisionHandlers(
        eqsys->GetHotTailCollisionHandler(),
        eqsys->GetRunawayCollisionHandler(),
        eqsys->GetREFluid()
    );
}


/**
 * Construct a SolverLinearlyImplicit object according to the
 * provided settings.
 *
 * s:    Settings object specifying how to construct the
 *       SolverLinearlyImplicit object.
 * u:    List of unknown quantities.
 * eqns: List of equations for the unknowns of the equation system.
 */
SolverLinearlyImplicit *SimulationGenerator::ConstructSolver_linearly_implicit(
    Settings* /*s*/, FVM::UnknownQuantityHandler *u,
    vector<UnknownQuantityEquation*> *eqns
) {
    return new SolverLinearlyImplicit(u, eqns);
}

/**
 * Construct a SolverSNES object according to the provided settings.
 *
 * s:    Settings object specifying how to construct the
 *       SolverLinearlyImplicit object.
 * u:    List of unknown quantities.
 * eqns: List of equations for the unknowns of the equation system.
 */
SolverSNES *SimulationGenerator::ConstructSolver_nonlinear_snes(
    Settings* /*s*/, FVM::UnknownQuantityHandler *u,
    vector<UnknownQuantityEquation*> *nontrivial_unknowns
) {
    return new SolverSNES(u, nontrivial_unknowns);
}

