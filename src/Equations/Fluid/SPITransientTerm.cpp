/**
 * Implementation of an Euler backward transient term multiplied
 * by an arbitrary grid-dependent weight function. 
 * Evaluation of weights must be implemented in derived classes. 
 */

#include <iostream>
#include "FVM/config.h"
#include "FVM/Matrix.hpp"
#include "DREAM/Equations/Fluid/SPITransientTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
SPITransientTerm::SPITransientTerm(FVM::Grid *grid, const len_t unknownId, len_t nShard, real_t scaleFactor)
        : FVM::DiagonalTerm(grid), unknownId(unknownId), nShard(nShard), scaleFactor(scaleFactor){}

/**
 * Rebuild the transient term.
 *
 * dt: Length of next time step to take.
 */
void SPITransientTerm::Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler *uqty) {
    this->dt = dt;
    this->xn = uqty->GetUnknownDataPrevious(this->unknownId);

    if(!hasBeenInitialized){
        InitializeWeights();
        hasBeenInitialized = true;
    }
    if(TermDependsOnUnknowns()) 
        SetWeights();
}

void SPITransientTerm::SetWeights() {
    for (len_t ip = 0; ip < nShard; ip++){
        weights[ip] = scaleFactor;
    }
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
void SPITransientTerm::SetMatrixElements(FVM::Matrix *mat, real_t *rhs) {
    for (len_t i = 0; i < nShard; i++)
        mat->SetElement(i, i, weights[i]/this->dt, ADD_VALUES);

    if (rhs != nullptr)
        for (len_t i = 0; i < nShard; i++)
            rhs[i] -= weights[i]*this->xn[i] / this->dt;
}

/**
 * Set the elements in the function vector 'F' (i.e.
 * evaluate this term).
 *
 * vec: Vector containing value of 'F' on return.
 * x:   Previous solution (unused).
 */
void SPITransientTerm::SetVectorElements(real_t *vec, const real_t *xnp1) {

    for (len_t i = 0; i < nShard; i++)
        vec[i] += weights[i]*(xnp1[i] - xn[i]) / this->dt;
}

