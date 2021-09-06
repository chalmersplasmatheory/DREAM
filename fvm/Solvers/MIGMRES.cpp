/**
 * Implementation of matrix invertor utilizing an iterative
 * Generalized Minimal Residual (GMRES) method with an LU preconditioner.
 */

#include <petscksp.h>
#include <petscvec.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Solvers/MIGMRES.hpp"

using namespace DREAM::FVM;
using namespace std;

/**
 * Constructor.
 *
 * n: Number of elements in solution vector.
 */
MIGMRES::MIGMRES(
    const len_t n, vector<len_t>& nontrivial_unknowns,
    UnknownQuantityHandler *unknowns
) {
    KSPCreate(PETSC_COMM_WORLD, &this->ksp);
    this->xn = n;

    // Construct a vector indicating to PETSc how to
    // partition the matrix into blocks (we divide the
    // matrix into one block per unknown quantity in
    // the system)
    this->nBlocks = nontrivial_unknowns.size();
    this->blocks = new PetscInt[this->nBlocks];
    // Each element of 'blocks' contains the size of
    // the corresponding matrix block
    for (len_t i = 0; i < this->nBlocks; i++) {
        this->blocks[i] =
            unknowns->GetUnknown(nontrivial_unknowns[i])
                ->NumberOfElements();
    }

    // Configure the preconditioner
    PC pc;
    KSPGetPC(ksp, &pc);

    PCSetType(pc, PCBJACOBI);
    PCBJacobiSetTotalBlocks(pc, this->nBlocks, this->blocks);
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

