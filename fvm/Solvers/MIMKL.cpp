/**
 * Implementation of matrix invertor based on LU
 * factorization (i.e. direct inversion).
 */

#include <omp.h>
#include <petscmat.h>
#include <petscvec.h>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Solvers/MIMKL.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MIMKL::MIMKL(const len_t n, bool verbose) {
    this->verbose = verbose;

    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;
}

/**
 * Destructor.
 */
MIMKL::~MIMKL() {
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
#ifdef PETSC_HAVE_MKL_PARDISO
void MIMKL::Invert(Matrix *A, Vec *b, Vec *x) {
    Mat F;
    PC pc;

    KSPSetOperators(this->ksp, A->mat(), A->mat());
    
    // Set direct LU factorization
    KSPGetPC(this->ksp, &pc);
    PCSetType(pc, PCLU);
    PCFactorSetMatSolverType(pc, MATSOLVERMKL_PARDISO);
    PCFactorSetUpMatSolverType(pc);
    PCFactorGetMatrix(pc, &F);
    KSPSetType(this->ksp, KSPPREONLY);

    // Force PARDISO to only use one CPU thread
    // (otherwise it will crash with an "insufficient memory"
    // error)
    MatMkl_PardisoSetCntl(F, 65, 1);

    // Maximum iterations for refinement
    MatMkl_PardisoSetCntl(F, 8, 10);

    if (this->verbose)
        MatMkl_PardisoSetCntl(F, 68, 1);

    // Solve
    this->errorcode = KSPSolve(this->ksp, *b, *x);

    if (this->errorcode != 0)
        PCPostSolve(pc, this->ksp);
#else
void MIMKL::Invert(Matrix*, Vec*, Vec*) {
#endif
}

