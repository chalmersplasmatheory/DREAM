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

    s->DefineSetting(MODULENAME "/backupsolver", "Type of backup linear solver to use if the main linear solver fails", (int_t)OptionConstants::LINEAR_SOLVER_NONE);
    s->DefineSetting(MODULENAME "/linsolv", "Type of linear solver to use", (int_t)OptionConstants::LINEAR_SOLVER_LU);
    s->DefineSetting(MODULENAME "/maxiter", "Maximum number of nonlinear iterations allowed", (int_t)100);
    s->DefineSetting(MODULENAME "/reltol", "Relative tolerance for nonlinear solver", (real_t)1e-6);
    s->DefineSetting(MODULENAME "/verbose", "If true, generates extra output during nonlinear solve", (bool)false);

    DefineToleranceSettings(MODULENAME, s);
    DefinePreconditionerSettings(s);

    // Debug settings
    s->DefineSetting(MODULENAME "/debug/printmatrixinfo", "Print detailed information about the PETSc matrix", (bool)false);
    s->DefineSetting(MODULENAME "/debug/printjacobianinfo", "Print detailed information about the jacobian PETSc matrix", (bool)false);
    s->DefineSetting(MODULENAME "/debug/savejacobian", "If true, saves the jacobian matrix in the specified iteration(s)", (bool)false);
    s->DefineSetting(MODULENAME "/debug/savesolution", "Saves the solution in the specified iteration, i.e. x = (J^-1)F", (bool)false);
    s->DefineSetting(MODULENAME "/debug/savematrix", "If true, saves the linear operator matrix in the specified time step(s)", (bool)false);
    s->DefineSetting(MODULENAME "/debug/savenumericaljacobian", "If true, evaluates the jacobian numerically and saves it for the specified iteration(s)", (bool)false);
    s->DefineSetting(MODULENAME "/debug/saverhs", "If true, saves the RHS vector in the specified iteration(s)", (bool)false);
    s->DefineSetting(MODULENAME "/debug/saveresidual", "If true, saves the residual vector in the specified iteration(s)", (bool)false);
    s->DefineSetting(MODULENAME "/debug/savesystem", "If true, saves the full equation system in the most recent iteration/time step", (bool)false);
    s->DefineSetting(MODULENAME "/debug/rescaled", "If true, saves the rescaled version of the jacobian/solution/residual.", (bool)false);
    s->DefineSetting(MODULENAME "/debug/timestep", "Index of time step to save debug info for. If '0', saves debug info for all time steps and iterations", (int_t)0);
    s->DefineSetting(MODULENAME "/debug/iteration", "Index of iteration to save debug info for.", (int_t)1);
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
            solver = ConstructSolver_linearly_implicit(s, u, eqns, eqsys);
            break;

		case OptionConstants::SOLVER_TYPE_NONLINEAR:
			solver = ConstructSolver_nonlinear(s, u, eqns, eqsys);
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

    solver->SetSPIHandler(eqsys->GetSPIHandler());

    solver->SetIonHandler(eqsys->GetIonHandler());

    solver->SetConvergenceChecker(LoadToleranceSettings(
        MODULENAME, s, eqns, u, solver->GetNonTrivials()
    ));
	solver->SetExternalIteratorConvergenceChecker(LoadToleranceSettings(
        MODULENAME, s, eqns, u, *eqsys->GetExternallyIteratedUnknowns()
	));

    solver->SetPreconditioner(LoadPreconditionerSettings(
        s, u, solver->GetNonTrivials()
    ));

	solver->SetExternalIteratorMaxIterations(s->GetInteger(MODULENAME "/maxiter"));
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
    vector<UnknownQuantityEquation*> *eqns,
    EquationSystem *eqsys
) {
    enum OptionConstants::linear_solver linsolv =
        (enum OptionConstants::linear_solver)s->GetInteger(MODULENAME "/linsolv");
    
    bool verbose    = s->GetBool(MODULENAME "/verbose");
    bool printdebug = s->GetBool(MODULENAME "/debug/printmatrixinfo");
    bool savematrix = s->GetBool(MODULENAME "/debug/savematrix");
    bool saverhs    = s->GetBool(MODULENAME "/debug/saverhs");
    int_t timestep  = s->GetInteger(MODULENAME "/debug/timestep");
    bool savesystem = s->GetBool(MODULENAME "/debug/savesystem");

    auto sli = new SolverLinearlyImplicit(u, eqns, eqsys, verbose, linsolv);
    sli->SetDebugMode(printdebug, savematrix, saverhs, timestep, savesystem);

    return sli;
}

/**
 * Construct a SolverNonLinear object according to the provided
 * settings.
 */
SolverNonLinear *SimulationGenerator::ConstructSolver_nonlinear(
	Settings *s, FVM::UnknownQuantityHandler *u,
	vector<UnknownQuantityEquation*> *eqns,
    EquationSystem *eqsys
) {
    enum OptionConstants::linear_solver
        backups = (enum OptionConstants::linear_solver)s->GetInteger(MODULENAME "/backupsolver"),
        linsolv = (enum OptionConstants::linear_solver)s->GetInteger(MODULENAME "/linsolv");

    int_t maxiter     = s->GetInteger(MODULENAME "/maxiter");
    real_t reltol     = s->GetReal(MODULENAME "/reltol");
    bool verbose      = s->GetBool(MODULENAME "/verbose");
    bool savejacobian = s->GetBool(MODULENAME "/debug/savejacobian");
    bool savesolution = s->GetBool(MODULENAME "/debug/savesolution");
    bool savenumjac   = s->GetBool(MODULENAME "/debug/savenumericaljacobian");
    bool saveresidual = s->GetBool(MODULENAME "/debug/saveresidual");
    bool printdebug   = s->GetBool(MODULENAME "/debug/printjacobianinfo");
    bool rescaled     = s->GetBool(MODULENAME "/debug/rescaled");
    int_t timestep    = s->GetInteger(MODULENAME "/debug/timestep");
    int_t iteration   = s->GetInteger(MODULENAME "/debug/iteration");
    bool savesystem   = s->GetBool(MODULENAME "/debug/savesystem");

    auto snl = new SolverNonLinear(u, eqns, eqsys, linsolv, backups, maxiter, reltol, verbose);
    snl->SetDebugMode(printdebug, savesolution, savejacobian, saveresidual, savenumjac, timestep, iteration, savesystem, rescaled);

    return snl;
}

