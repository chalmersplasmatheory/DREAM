/**
 * Implementation of matrix invertor based on LU
 * factorization (i.e. direct inversion).
 */

#include <petscmat.h>
#include <petscvec.h>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Solvers/MISuperLU.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MISuperLU::MISuperLU(const len_t n) {
    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;
}

/**
 * Destructor.
 */
MISuperLU::~MISuperLU() {
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
void MISuperLU::Invert(Matrix *A, Vec *b, Vec *x) {
    PC pc;

    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Set direct LU factorization
    KSPGetPC(this->ksp, &pc);
    PCSetType(pc, PCLU);
    PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU);
    KSPSetType(this->ksp, KSPPREONLY);

    // Solve
    KSPSolve(this->ksp, *b, *x);
}

