/**
 * Implementation of matrix invertor based on LU
 * factorization (i.e. direct inversion).
 */

#include <petscvec.h>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Solvers/MIMUMPS.hpp"

using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MIMUMPS::MIMUMPS(const len_t n) {
    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;
}

/**
 * Destructor.
 */
MIMUMPS::~MIMUMPS() {
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
void MIMUMPS::Invert(Matrix *A, Vec *b, Vec *x) {
    PC pc;
    Mat F;

    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Set direct LU factorization
    KSPGetPC(this->ksp, &pc);
    PCSetType(pc, PCLU);
    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
    //PCFactorSetUpMatSolverType(pc);
    KSPSetType(this->ksp, KSPPREONLY);

    // Solve
    KSPSolve(this->ksp, *b, *x);

    PCFactorGetMatrix(pc, &F);

    PetscInt info1, info2;
    MatMumpsGetInfo(F, 1, &info1);
    MatMumpsGetInfo(F, 2, &info2);

    if (info1 != 0) {
        printf(":: MUMPS INFO(1) = %d\n", info1);
        printf(":: MUMPS INFO(2) = %d\n", info2);
    }
}

