/**
 * Implementation of the 'SNES_set_jacobian()' and 'SNES_set_function()'
 * routines, which builds the jacobian matrix and function vector for
 * the SNES solver.
 */

#include <cstdio>
#include <petsc.h>
#include "DREAM/IO.hpp"
#include "DREAM/Solver/SolverSNES.hpp"


using namespace DREAM;


/**
 * Custom convergence test for the non-linear solver.
 *
 * snes:   SNES context.
 * it:     Current iteration (0 is the first is before any Newton step).
 * xnorm:  2-norm of current iterate.
 * gnorm:  2-norm of current step.
 * f:      2-norm of function.
 * reason: On return, contains the 'converged reason' which indicates
 *         to SNES whether the solution is converged, diverged, or needs
 *         further iterations.
 * cctx:   Convergence context (a pointer to the parent 'SolverSNES' object).
 */
PetscErrorCode DREAM::SNES_convergence_test(
    SNES snes, PetscInt it, PetscReal /*xnorm*/,
    PetscReal /*gnorm*/, PetscReal /*f*/,
    SNESConvergedReason *reason, void *cctx
) {
    SolverSNES *solver = (SolverSNES*)cctx;

    // Always require at least one iteration
    if (it == 0) {
        *reason = SNES_CONVERGED_ITERATING;
        return 0;
    // Guard for excessive number of iterations
    } else if (it >= solver->MaxIter()) {
        *reason = SNES_DIVERGED_MAX_IT;
        return 0;
    }

    // Get dx from SNES
    Vec petsc_x, petsc_dx;
    SNESGetSolution(snes, &petsc_x);
    SNESGetSolutionUpdate(snes, &petsc_dx);

    // Obtain 2-norms separately for each unknown quantity
    real_t *x, *dx;
    VecGetArray(petsc_x, &x);
    VecGetArray(petsc_dx, &dx);

    // Get vectors to store norms in
    real_t *x_2norm  = solver->Get_x_2norm();
    real_t *dx_2norm = solver->Get_dx_2norm();

    solver->CalculateNonTrivial2Norm(x, x_2norm);
    solver->CalculateNonTrivial2Norm(dx, dx_2norm);

    VecRestoreArray(solver->GetPETScSolution(), &x);
    VecRestoreArray(petsc_dx, &dx);

    // Iterate over norms and ensure that all are small
    const len_t N = solver->GetNUnknowns();
    bool converged = true;
    const real_t stol = solver->RelTol();

    if (solver->Verbose())
        DREAM::IO::PrintInfo("ITERATION %d", it);

    for (len_t i = 0; i < N; i++) {
        bool conv = (dx_2norm[i] < stol*x_2norm[i]);

        if (solver->Verbose()) {
#ifdef COLOR_TERMINAL
            if (conv)
                DREAM::IO::PrintInfo(
                    "   \x1B[32m%10s  |x| = %e, |dx| = %e\x1B[0m",
                    solver->GetNonTrivialName(i).c_str(),
                    x_2norm[i], dx_2norm[i]
                );
             else
                DREAM::IO::PrintInfo(
                    "   \x1B[1;31m%10s  |x| = %e, |dx| = %e\x1B[0m",
                    solver->GetNonTrivialName(i).c_str(),
                    x_2norm[i], dx_2norm[i]
                );
#else
            DREAM::IO::PrintInfo(
                "   %10s  |x| = %e, |dx| = %e",
                solver->GetNonTrivialName(i).c_str(),
                x_2norm[i], dx_2norm[i]
            );
#endif
        }

        converged = converged && conv;
    }

    // Determine if solver has converged
    if (converged)
        *reason = SNES_CONVERGED_SNORM_RELATIVE;
    else
        *reason = SNES_CONVERGED_ITERATING;

    return 0;
}

/**
 * Build the jacobian matrix for the SNES solver.
 *
 * snes: SNES context.
 * x:    Vector to evaluate jacobian with.
 * Amat: (Approximate) jacobian matrix (same as 'solver->jacobian')
 * Pmat: The matrix to be used in constructing the preconditioner,
 *       usually the same as Amat.
 * ctx:  A 'SolverSNES' object.
 */
PetscErrorCode DREAM::SNES_set_jacobian(
    SNES snes, Vec /*x*/, Mat /*Amat*/, Mat /*Pmat*/, void *ctx
) {
    SolverSNES *solver = (SolverSNES*)ctx;

    PetscInt iter;
    SNESGetIterationNumber(snes, &iter);

    // Rebuild equation terms (if necessary)
    SNES_update_system(iter, solver);
    solver->BuildJacobian(solver->CurrentTime(), solver->CurrentTimeStep(), solver->GetJacobian());

    return 0;
}

/**
 * Build the function vector for the SNES solver.
 * 
 * snes: SNES context.
 * x:    State at which to evaluate residual.
 * f:    Vector to store function value in.
 * ctx:  A 'SolverSNES' object.
 */
PetscErrorCode DREAM::SNES_set_function(SNES snes, Vec /*x*/, Vec f, void *ctx) {
    SolverSNES *solver = (SolverSNES*)ctx;

    PetscInt iter;
    SNESGetIterationNumber(snes, &iter);

    // Rebuild equation terms (if necessary)
    SNES_update_system(iter, solver);

    real_t *fvec;
    VecGetArray(f, &fvec);

    solver->BuildVector(solver->CurrentTime(), solver->CurrentTimeStep(), fvec, solver->GetJacobian());

    VecRestoreArray(f, &fvec);

    return 0;
}

/**
 * Rebuild all equation terms (ahead of building the function
 * vector and jacobian).
 *
 * iter:   Current iteration number.
 * solver: DREAM Solver object running the SNES solve.
 */
void DREAM::SNES_update_system(const int_t /*iter*/, SolverSNES *solver) {
    /*if (solver->GetLastRebuild() < iter) {
        solver->SetLastRebuild(iter);
        solver->RebuildTerms(solver->CurrentTime(), solver->CurrentTimeStep());
    }*/
    solver->RebuildTerms(solver->CurrentTime(), solver->CurrentTimeStep());
}

/**
 * Function which is called after each Newton iteration when an
 * updated solution is available. This function inserts the resulting
 * data into the various 'QuantityData' objects.
 *
 * snes: SNES context.
 * its:  Iteration number.
 * norm: 2-norm function value.
 * mctx: Monitoring context (the DREAM 'SolverSNES' object here).
 */
PetscErrorCode DREAM::SNES_solution_obtained(
    SNES /*snes*/, PetscInt its, PetscReal /*norm*/, void *mctx
) {
    SolverSNES *solver = (SolverSNES*)mctx;
    solver->StoreSolution(its);

    return 0;
}

