/**
 * Implementation of an Euler backward transient term multiplied
 * by an arbitrary grid-dependent weight function. 
 * Evaluation of weights must be implemented in derived classes. 
 */

#include <iostream>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "FVM/Equation/WeightedTransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM::FVM;


/**
 * Constructor.
 */
WeightedTransientTerm::WeightedTransientTerm(Grid *grid, const len_t unknownId)
        : EquationTerm(grid), unknownId(unknownId){ 
    //weights = new real_t[grid->GetNCells()];
    //GridRebuilt();
    //if(!TermDependsOnUnknowns())
    //    SetWeights();
}

/**
 * Destructor
 */
WeightedTransientTerm::~WeightedTransientTerm(){
    this->DeallocateMemory();
    if(weights!=nullptr)
        delete [] weights;
    }


/**
 * Called if the grid is rebuilt; reallocates and rebuilds quantities.
 */
bool WeightedTransientTerm::GridRebuilt(){
    this->AllocateMemory();
    if(weights!=nullptr)
        delete [] weights;
    weights = new real_t[grid->GetNCells()];
    if(!TermDependsOnUnknowns())
        SetWeights();
    return true;
}

/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 */
void WeightedTransientTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, Matrix *mat, const real_t* /*x*/
) {
    if (uqtyId == derivId)
        this->SetMatrixElements(mat, nullptr);
}

/**
 * Set the matrix elements corresponding to this
 * transient term.
 *
 * This term assumes that the linearized matrix equation
 * to solve is of the form
 *
 *   df/dt + Mf = -S
 *
 * where M is a linear matrix operator represented by 'mat',
 * and 'S' is a source term stored in 'rhs'.
 *
 * mat: Matrix to set elements of.
 * rhs: Equation RHS.
 */
void WeightedTransientTerm::SetMatrixElements(Matrix *mat, real_t *rhs) {
    const len_t N = grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]/this->dt, ADD_VALUES);

    if (rhs != nullptr)
        for (len_t i = 0; i < N; i++)
            rhs[i] += weights[i]*this->xn[i] / this->dt;
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution (unused).
 */
void WeightedTransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {
    const len_t N = grid->GetNCells();

    for (len_t i = 0; i < N; i++)
        vec[i] += weights[i]*(xnp1[i] - xn[i]) / this->dt;
}

