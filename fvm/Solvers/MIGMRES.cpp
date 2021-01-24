/**
 * Implementation of matrix invertor utilizing an iterative
 * Krylov Subspace (KSP) method.
 */

#include <petscvec.h>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Solvers/MIKSP.hpp"

using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MIKSP::MIKSP(const len_t n) {
    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;
}

/**
 * Destructor.
 */
MIKSP::~MIKSP() {
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
void MIKSP::Invert(Matrix *A, Vec *b, Vec *x) {
    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Solve
    KSPSetType(this->ksp, KSPGMRES);
    KSPSolve(this->ksp, *b, *x);
}

/*'
 * Set the function to use for checking if the
 * solution is converged.
 */
void MIKSP::SetConvergenceTest(
    PetscErrorCode (*converge)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*),
    void *context
) {
    KSPSetConvergenceTest(this->ksp, converge, context, nullptr);
}

