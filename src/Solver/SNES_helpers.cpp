/**
 * Implementation of the 'SNES_set_jacobian()' and 'SNES_set_function()'
 * routines, which builds the jacobian matrix and function vector for
 * the SNES solver.
 */

#include <cstdio>
#include <petsc.h>
#include "DREAM/Solver/SolverSNES.hpp"


using namespace DREAM;


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
PetscErrorCode DREAM::SNES_set_jacobian(SNES /*snes*/, Vec /*x*/, Mat /*Amat*/, Mat /*Pmat*/, void *ctx) {
    SolverSNES *solver = (SolverSNES*)ctx;
    printf("[SNES] Evaluate jacobian\n");

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
PetscErrorCode DREAM::SNES_set_function(SNES /*snes*/, Vec x, Vec f, void *ctx) {
    SolverSNES *solver = (SolverSNES*)ctx;
    printf("[SNES] Evaluate function\n");

    // Rebuild equation terms
    SNES_update_system(solver);

    /*PetscInt idx = 6;
    PetscScalar v;
    VecGetValues(x, 1, &idx, &v);
    printf("f0 = %e\n", v);*/

    real_t *fvec;
    VecGetArray(f, &fvec);
    solver->BuildVector(solver->CurrentTime(), solver->CurrentTimeStep(), fvec, solver->GetJacobian());

    //printf("f0 = %e\n", fvec[6]);
    
    VecRestoreArray(f, &fvec);

    return 0;
}

/**
 * Rebuild all equation terms (ahead of building the function
 * vector and jacobian).
 *
 * solver: DREAM Solver object running the SNES solve.
 */
void DREAM::SNES_update_system(SolverSNES *solver) {
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
PetscErrorCode DREAM::SNES_solution_obtained(SNES /*snes*/, PetscInt its, PetscReal norm, void *mctx) {
    printf("Iteration %d, |x| = %.12e\n", its, norm);

    // TODO
    SolverSNES *solver = (SolverSNES*)mctx;
    solver->StoreSolution(its);

    return 0;
}

