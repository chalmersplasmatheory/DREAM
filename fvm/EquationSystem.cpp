/**
 * Implementation of the DREAM equation system class.
 */

#include <iostream>
#include <petscmat.h>
#include "FVM/config.h"
#include "FVM/EquationSystem.hpp"

using namespace DREAM::FVM;
using namespace std;

/**
 * Constructor.
 */
EquationSystem::EquationSystem() { }

/**
 * Destructor.
 */
EquationSystem::~EquationSystem() {
    this->Destroy();
}

/**
 * Construct the matrix
 */
void EquationSystem::ConstructSystem() {
    // Determine matrix size
    PetscInt mSize = this->next_subindex;

    // Calculate number of non-zero elements in matrix
    PetscInt nnz = 0;
    for (vector<struct _subeq>::iterator it = this->subeqs.begin(); it != this->subeqs.end(); it++) {
        nnz += it->nnz;
    }

    this->Construct(mSize, mSize, nnz);
}

/**
 * Defines a new "sub-equation" to include in the matrix. The
 * sub-equation appears by getting its own (square) matrix block
 * in the large matrix, and can be accessed as a separate matrix.
 *
 * n:    Number of elements in solution vector.
 * nnz:  Number of non-zero elements in matrix block row
 *       representing this equation.
 */
len_t EquationSystem::CreateSubEquation(const PetscInt n, const PetscInt nnz) {
    // Define index set
    struct _subeq se;
    se.n      = n;
    se.nnz    = nnz;
    se.offset = this->next_subindex;

    this->subeqs.push_back(se);

    this->next_subindex += n;

    return (this->subeqs.size()-1);
}

/**
 * Sets which sub-equation to write into.
 *
 * subeq1: Index of equation to set matrix to (block row index of sub-matrix).
 * subeq2: Index of unknown to set matrix to (block row column of sub-matrix).
 */
void EquationSystem::SelectSubEquation(const PetscInt subeq1, const PetscInt subeq2) {
    this->SetOffset(this->subeqs.at(subeq1).offset, this->subeqs.at(subeq2).offset);
}

/**
 * Sets the matrix corresponding to the specified
 * sub-equation to zero.
 *
 * subeq1: Index of equation for which matrix should be zeroed (block row index of sub-matrix).
 * subeq2: Index of unknown for which matrix should be zeroed (block row column of sub-matrix).
 */
void EquationSystem::ZeroEquation(const PetscInt subeq) {
    IS is;
    ISCreateStride(PETSC_COMM_WORLD, this->subeqs.at(subeq).n, this->subeqs.at(subeq).offset, 1, &is);

    MatZeroRowsColumnsIS(this->petsc_mat, is, 0, nullptr, nullptr);

    ISDestroy(&is);
}

