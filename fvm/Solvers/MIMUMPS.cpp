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
    #define CHKERR(e) do { if (ierr != 0) throw FVMException("MUMPS error: %d", ierr); } while (false)
    PetscErrorCode ierr;
    PC pc;

    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Set direct LU factorization
    KSPGetPC(this->ksp, &pc);
    PCSetType(pc, PCLU);
    ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS); CHKERR(ierr);
    //ierr = PCFactorSetUpMatSolverType(pc); CHKERR(ierr);      // Needed for setting control parameters
    KSPSetType(this->ksp, KSPPREONLY);

    // Set MUMPS control parameters
    //Mat F;
    //ierr = PCFactorGetMatrix(pc, &F); CHKERR(ierr);
    //MatMumpsSetIcntl(F, ICNTL_OPENMP_THREADS, 2);
    //MatMumpsSetIcntl(F, ICNTL_PRINTING_LEVEL, 2);
    //MatMumpsSetIcntl(F, ICNTL_DETECT_NULL_PIVOT_ROWS, 1);
    //MatMumpsSetCntl(F, 1, 0.0);   // Disable pivoting

    // Solve
    ierr = KSPSolve(this->ksp, *b, *x); CHKERR(ierr);
}

