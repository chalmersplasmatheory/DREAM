/**
 * Construct a time stepper object.
 */

#include "DREAM/IO.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/Solver/SolverLinearlyImplicit.hpp"
#include "DREAM/Solver/SolverNonLinear.hpp"
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
    s->DefineSetting(MODULENAME "/type", "Equation system solver type", (int_t)OptionConstants::SOLVER_TYPE_NONLINEAR);

    s->DefineSetting(MODULENAME "/linsolv", "Type of linear solver to use", (int_t)OptionConstants::LINEAR_SOLVER_LU);
    s->DefineSetting(MODULENAME "/maxiter", "Maximum number of nonlinear iterations allowed", (int_t)100);
    s->DefineSetting(MODULENAME "/reltol", "Relative tolerance for nonlinear solver", (real_t)1e-6);
    s->DefineSetting(MODULENAME "/verbose", "If true, generates extra output during nonlinear solve", (bool)false);

    DefineToleranceSettings(MODULENAME, s);
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

		case OptionConstants::SOLVER_TYPE_NONLINEAR:
			solver = ConstructSolver_nonlinear(s, u, eqns);
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

    ConvergenceChecker *cc = LoadToleranceSettings(
        MODULENAME, s, u, solver->GetNonTrivials()
    );

    solver->SetConvergenceChecker(cc);
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
    Settings *s, FVM::UnknownQuantityHandler *u,
    vector<UnknownQuantityEquation*> *eqns
) {
    enum OptionConstants::linear_solver linsolv =
        (enum OptionConstants::linear_solver)s->GetInteger(MODULENAME "/linsolv");

    return new SolverLinearlyImplicit(u, eqns, linsolv);
}

/**
 * Construct a SolverNonLinear object according to the provided
 * settings.
 */
SolverNonLinear *SimulationGenerator::ConstructSolver_nonlinear(
	Settings *s, FVM::UnknownQuantityHandler *u,
	vector<UnknownQuantityEquation*> *eqns
) {
    enum OptionConstants::linear_solver linsolv =
        (enum OptionConstants::linear_solver)s->GetInteger(MODULENAME "/linsolv");
    int_t maxiter = s->GetInteger(MODULENAME "/maxiter");
    real_t reltol = s->GetReal(MODULENAME "/reltol");
    bool verbose  = s->GetBool(MODULENAME "/verbose");

    return new SolverNonLinear(u, eqns, linsolv, maxiter, reltol, verbose);
}

