/**
 * Implementation of matrix invertor utilizing an iterative
 * Generalized Minimal Residual (GMRES) method with an LU preconditioner.
 */

#include <petscvec.h>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Solvers/MIGMRES.hpp"

using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MIGMRES::MIGMRES(const len_t n) {
    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;
}

/**
 * Destructor.
 */
MIGMRES::~MIGMRES() {
    KSPDestroy(&this->ksp);
    VecDestroy(&this->x);

    delete [] this->x_data;
}

/**
 * Solves the linear equation system represented by
 *
 *   Ax = b
 *
 * where A is a matrix, and b and x are vectors. A
 * pointer is returned to the solution, x.
 *
 * A: Matrix of size m-by-n representing the linear system.
 * b: Right-hand-side vector containing n elements.
 * x: Solution vector. Contains solution on return. Must be
 *    of size n at least.
 */
void MIGMRES::Invert(Matrix *A, Vec *b, Vec *x) {
    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Solve
    KSPSetType(this->ksp, KSPGMRES);
    KSPSolve(this->ksp, *b, *x);
}

/*'
 * Set the function to use for checking if the
 * solution is converged.
 */
void MIGMRES::SetConvergenceTest(
    PetscErrorCode (*converge)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*),
    void *context
) {
    KSPSetConvergenceTest(this->ksp, converge, context, nullptr);
}

